#include <Arduino.h>
#include <map>
#include <string>
#include <Wire.h> 
#include <LiquidCrystal_I2C.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>

//web server stuff
#include <SPIFFS.h>
#include <WiFi.h>
#include <ESPAsyncWebServer.h>
#include <AsyncTCP.h>
#include <ArduinoJson.h>

const char* ssid = "ESP32_BWRSim";
const char* password = "";

AsyncWebServer server(80);
AsyncWebSocket ws("/ws");

LiquidCrystal_I2C lcd(0x27,20,4);
Adafruit_SSD1306 display(128, 64, &Wire, -1);
//constants
struct CrossSections {
  byte U235_capture = 99;
  int U235_fission = 583;
  byte U238_capture = 2;
  byte H2_scatter = 4;
  double C12_capture = 0.002;
  int B10_capture = 3840;
};
const CrossSections crossSections;
struct Direction {
  short x;
  short y;
};
const Direction directions[] = {
  {6,0},
  {-6,0},
  {0,6},
  {0,-6}
};
const double BARN = 1e-24;
const double N235 = 2.53e22;
const double N238 = 2.51e22;
const double NB10 = 5e21;
const double FullPowerFlux = 2.5e12;
//not constants
struct ReactorInventory {
  int waterMass = 75000;
};
ReactorInventory reactorInventory;
struct Model {
  double reactorWaterTemperature = 20;
  double RRCAPercent = 0.05;
  double RRCBPercent = 0.05;
};
Model model;
const char* rodPositions[] = {
                        "15-03","21-03","27-03",
                "09-09","15-09","21-09","27-09","33-09",
        "03-15","09-15","15-15","21-15","27-15","33-15","39-15",
        "03-21","09-21","15-21","21-21","27-21","33-21","39-21",
        "03-27","09-27","15-27","21-27","27-27","33-27","39-27",
                "09-33","15-33","21-33","27-33","33-33",
                        "15-39","21-39","27-39"
};
struct Fuel { //return type garbage for fuel_get
  double kStep;
  double kEff;
};
enum RodMotion {
  IDLE,
  INSERT,
  WITHDRAW,
  SETTLE
};
struct RodState {
  double n = 100; //local neutron flux
  byte p; //control rod position (0 to 48)
  RodMotion state = IDLE;

  RodState(double neutronFlux = 100, byte pos = 0, RodMotion s = IDLE)  : n(neutronFlux), p(pos), state(s) {}
};
std::map<std::string, RodState> rods;
Fuel fuel_get(double dt, double waterMass, double controlDepth, double neutronFlux, double temperatureFuel, double coreFlow){
  double SigmaF235 = N235 * crossSections.U235_fission * BARN;
  double SigmaC235 = N235 * crossSections.U235_capture * BARN;
  double SigmaC238 = N238 * crossSections.U238_capture * BARN;

  double SigmaAFuel = SigmaF235 + SigmaC235 + SigmaC238;

  double v = 2.43;
  double eta = (v * SigmaF235) / SigmaAFuel;

  double U = SigmaAFuel;
  double M = crossSections.C12_capture * BARN * N235;

  double insertion = 1 - controlDepth;
  double CR = (crossSections.B10_capture * BARN * NB10) * insertion;

  double f = U / (U + M + CR);

  double voids = neutronFlux / FullPowerFlux;
  double totalCoreFlow = constrain(abs(coreFlow-1),0,1);
  voids = constrain(totalCoreFlow*voids,0,0.7);

  f = f*(1-voids);

  double p = 0.9;
  double eps = 1.03;
  double Pfn1 = 0.97;
  double Ptn1 = 0.97;

  double alpha = 0.0001;
  byte Tref = 20;

  double kEff = eta * f * p * eps * Pfn1 * Ptn1;
  kEff *= (1 - alpha * (temperatureFuel - Tref));
  double kStep = pow(kEff,dt * 5);

  return {kStep, kEff};
};

void reactor_physics(double dt, std::map<std::string, RodState>& rods){
  double totalCoreFlow = model.RRCAPercent + model.RRCBPercent;
  double waterMass = reactorInventory.waterMass;

  std::map<std::string, double> produced_n; //name = neutron flux
  std::map<std::string, double> k_factors; //name = kStep

  double avg_power = 0;
  double old_temp = model.reactorWaterTemperature;
  double new_temp = model.reactorWaterTemperature;
  for (auto& [name, info]: rods){
    double NeutronFlux = info.n; if (NeutronFlux < 100){NeutronFlux = 100;} // "local NeutronFlux = math.max(info.n, 100) because C++ says FUCK YOU"
    Fuel mykEffArgs = fuel_get(dt, waterMass, double(info.p)/48, NeutronFlux, model.reactorWaterTemperature, totalCoreFlow);

    k_factors[name] = mykEffArgs.kStep;
    produced_n[name] = info.n * mykEffArgs.kStep; if (produced_n[name] < 10){produced_n[name] = 10;} //you get the point
  
    double energy = info.n/FullPowerFlux;
    avg_power += energy;
    energy = energy*760; //*delta # in MWt

    double calories = (energy*1000000)/rods.size();
    double HeatC = calories/1000;

    double TempNow = (HeatC/waterMass)*dt;
    new_temp += TempNow;
  }

  avg_power /= rods.size();

  std::map<std::string, double> diffusion_deltas;

  for (auto& [name, info]: rods){
    byte x = atoi(name.substr(0,2).c_str());
    byte y = atoi(name.substr(3,2).c_str());
    double my_amount = produced_n[name];

    for (const Direction& dir : directions){
      byte nx = x + dir.x; std::string nxs = "";
      byte ny = y + dir.y; std::string nys = "";
      if (std::to_string(nx).length() == 1){
        nxs = "0" + std::to_string(nx);
      }else{
        nxs = std::to_string(nx);
      }
      if (std::to_string(ny).length() == 1){
        nys = "0" + std::to_string(ny);
      }else{
        nys = std::to_string(ny);
      }
      std::string neighbor_name = nxs + "-" + nys;

      if (produced_n.find(neighbor_name) != produced_n.end()){
        double neighbor_amount = produced_n[neighbor_name];
        double transport_rate = k_factors[name];
        double transported = (my_amount - neighbor_amount) / 10 * transport_rate * dt;

        diffusion_deltas[name] -= transported;
        diffusion_deltas[neighbor_name] += transported;
      }else{
        double transport_rate = k_factors[name];

        double lost_flux = my_amount * 0.94;
        double transported = (my_amount - lost_flux)/10 * transport_rate;

        diffusion_deltas[name] -= transported;
      }
    }
  }
  for (auto& [name,val] : produced_n){
    double diff = diffusion_deltas[name] / 2;
    double final = val + diff; if (final < 10){final = 10;}

    rods[name].n = final;
  }

  double time_since_sd = 0;
  double power_before_sd = 0;

  if (avg_power < 0.02){
    if (time_since_sd == 0){
      power_before_sd = avg_power;
    }
    time_since_sd += dt;
  }else{
    power_before_sd = 0;
    time_since_sd = 0;
  }
  double pw = power_before_sd;
  double t_0 = 100*86400;
  double t = time_since_sd+1;

  double decay = 0.0622 * pw * (pow(t,-0.2)-pow((t_0 + t),-0.2));
  decay *= 1.3;

  double heat_generated = decay*760;
  double calories = heat_generated*1000000;

  double HeatC = calories/1000;

  double TempNow = (HeatC/waterMass * 4.18)*dt;
  new_temp += TempNow;

  const byte boilTemp = 110; //so much for the core being slightly pressurized smh
  if (model.reactorWaterTemperature > boilTemp){
    double excess = model.reactorWaterTemperature - boilTemp;

    double steam_generated = excess * 1;
    model.reactorWaterTemperature -= steam_generated * 0.01;
  }

  model.reactorWaterTemperature += new_temp - old_temp;
}

void setup() {
  Serial.begin(9600);
  SPIFFS.begin(true);
  WiFi.softAP(ssid, password);
  Serial.println(WiFi.softAPIP());
  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send(SPIFFS, "/web/index.html", "text/html");
  });

  server.addHandler(&ws);
  server.begin();

  lcd.init();
  lcd.backlight();
  lcd.clear();  
  for (const char* rod : rodPositions){
    Serial.println(rod);
    rods[std::string(rod)] = {100, 0, IDLE};
  }
  server.serveStatic("/", SPIFFS, "/web/");
  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send(SPIFFS, "/web/index.html", "text/html");
  });
}

String header;
unsigned long lastSent = 0;

void sendReactorData(double reactorPower){
  StaticJsonDocument<2048> doc;

  doc["APRM"] = reactorPower;
  doc["RRCA"] = model.RRCAPercent;
  doc["RRCB"] = model.RRCBPercent;

  JsonObject rodsObj = doc.createNestedObject("rods");

  for (auto& [name,info] : rods) {
    JsonObject rod = rodsObj.createNestedObject(name.c_str());
    rod["p"] = info.p;
  }

  String json;
  serializeJson(doc, json);

  ws.textAll(json);
}

double TimeSinceLastMovement = 0;
bool scramming = false;
bool doneScramming = true;
double TimeSinceLastScramMovement = 0;
bool lastMovementWasWithdraw = false;

void CDRM(float dt){
  TimeSinceLastMovement += dt;

  if (scramming) return;
  int joystickValue = analogRead(34);

  if (joystickValue >= 0 && joystickValue <= 2000){
    if (TimeSinceLastMovement < 1) return;
    TimeSinceLastMovement = 0;
    lastMovementWasWithdraw = true;
    for (auto& [name, info] : rods){
      info.state = WITHDRAW;
      if (info.p < 48) info.p++;
    }
  }
  else if (joystickValue > 3500) {
    if (TimeSinceLastMovement < 1) return;
    TimeSinceLastMovement = 0;
    lastMovementWasWithdraw = false;
    for (auto& [name, info] : rods){
      info.state = INSERT;
      if (info.p > 0) info.p--;
    }
  }
  else {
    for (auto & [name, info] : rods){
      if (info.state == WITHDRAW || info.state == INSERT) {
        if (info.p % 2 != 0){
          info.state = SETTLE;
        }
        else {
          info.state = IDLE;
        }
      }
    }
  }
}

float settleTimer = 0;

void CDRMSETTLE(float dt){
  settleTimer += dt;

  if (settleTimer < 0.2) return;
  settleTimer = 0;

  for (auto& [name, info] : rods){
    if (info.state == SETTLE){
      if (info.p % 2 != 0){
        if (lastMovementWasWithdraw) info.p++;
        if (!lastMovementWasWithdraw) info.p--;
      }
      if (info.p % 2 == 0){
        info.state = IDLE;
      }
    }
  }
}

void SCRAM(float dt){
  TimeSinceLastScramMovement += dt;

  if (not scramming) return;
  if (TimeSinceLastScramMovement <= 0.1) return;

  TimeSinceLastScramMovement = 0;

  float avgRodPos = 0.00;
  for (auto& [name,info] : rods){
    avgRodPos += info.p;
  }
  avgRodPos /= rods.size();

  if (!doneScramming and avgRodPos > 0.01) {
    for (auto& [name,info] : rods){
      if (info.p > 0) info.p--;
    }
  }
  else{
    doneScramming = true;
    scramming = false;
  }
}

void loop() {
  float dt = 5.00/60.00;

  if (digitalRead(2) == HIGH && !scramming){
    scramming = true;
    doneScramming = false;
  }
  CDRM(dt);
  SCRAM(dt);
  CDRMSETTLE(dt);

  reactor_physics(dt,rods);
  double neutrons = 0;
  for (auto& [name,info] : rods){
    neutrons += info.n;
  }
  double avg_flux = neutrons/ (double) rods.size();
  double APRM = (avg_flux/FullPowerFlux)*100.00;
  if (millis() - lastSent >= 200) {
    lastSent = millis();
    sendReactorData(APRM);
  }
  Serial.println(String(neutrons) + " C/s | " + String(APRM) + " %NOM | " + String(model.reactorWaterTemperature) + " Deg C | " + rods["15-03"].p + " Notches");
  std::string rodpos;
  if (std::to_string(rods["15-03"].p).length() == 1){
    rodpos = std::string("0") + std::to_string(rods["15-03"].p);
  }else{
    rodpos = std::to_string(rods["15-03"].p);
  }
  lcd.setCursor(0,0);
  lcd.print("R_p:");
  lcd.print(rodpos.c_str());
  lcd.print(" N_RE:");
  lcd.print(APRM);
  lcd.print("%");
  lcd.setCursor(0,1);
  lcd.print("T_POskut:");
  lcd.print(String(model.reactorWaterTemperature));
  lcd.print(" C");
  sleep(dt);
}