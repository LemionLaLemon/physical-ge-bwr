# Physical BWR/5 Simulator (simplified for my sanity)
This project is a physical simulation of BWR/5 reactors, particularly ones with a Mk 2 containment like Columbia Generating Station or LaSalle County Nuclear Generating Station.

![LaSalle County Nuclear Generating Station](https://www.chicagotribune.com/wp-content/uploads/migration/2014/12/09/LWAGENEIPRGF5OAEISVTYSRXCE.jpg?w=1200&resize=1200,900)
![Columbia Generating Station](https://thumb.spokesman.com/fmgC7__1Ka2-X9CiQJ9qpVjCQ5g=/1200x800/smart/media.spokesman.com/photos/2011/05/04/2.jpg)

# Firmware
This project uses physics explained in the [R304B](https://www.nrc.gov/docs/ML1125/ML11258A296.pdf) document package, and code examples were taken from [OpenLWR](https://github.com/OpenLWR/)

## fuel_get()
This section is the fuel_get function that returns local kEff (reactivity multiplier) and kStep (speed) from the cross sections, control rod position, recirculation flow, water mass, and water temperature. This function was derived from the 6 term kEff equation found in R304B.

![K_eff equation](./readme%20images/keff.png)

where ε is the fast fission factor, L_f is the fast nonleakage factor, ρ is the resonance escape probability, L_th is the thermal nonleakage factor, f is the thermal utilization factor made up from sums of cross sections for absorption, and η is the reproduction factor made up from sums of the number of neutrons per fission.

f can be calculated using this equation

![f equation](./readme%20images/f.png)

and η can be calculated using this equation

![n equation](./readme%20images/n.png)

The kEff and kStep calculated in this function can then be used to calculate the local neutron population for each rod in the next iteration.

```cpp
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
  double totalCoreFlow = constrain(abs((coreFlow/100)-1),0,1);
  voids = constrain(totalCoreFlow*voids,0,0.7);

  f = f*(1-voids);

  double p = 0.9;
  double eps = 1.03;
  double Pfn1 = 0.97;
  double Ptn1 = 0.97;

  double alpha = 0.0001;
  double Tref = 20;

  double kEff = eta * f * p * eps * Pfn1 * Ptn1;
  kEff *= (1 - alpha * (temperatureFuel - Tref));
  double kStep = pow(kEff,dt * 5);

  return {kStep, kEff};
};
```

## reactor_physics()
This section uses the kEff and kStep from fuel_get() to calculate the amount of neutrons for each control rod for the next generation. This can now be simply calculated as just n_1 = n_0 * kEff

>If there are (n_0) neutrons in one generation, then there will be (n_0 * K_eff) neutrons present in the next generation.
\
>source: [R304B, page 12](https://www.nrc.gov/docs/ML1125/ML11258A296.pdf)

This section is only as long as it is due to the part that calculates neutron diffusion. When neutrons diffuse to neighboring rods as neutrons do.

![Neutron Diffusion Visual](./readme%20images/diffusion.png)

In the image above, the rod with the most neutrons (the red rod) will diffuse 1/4th of its neutrons to each of its neighbor. Then, each neighbor will diffuse 1/4 of their neutrons to their neighbor, and so on. This is also why the center of the core is usually the section with the highest local power, since there's more neighbors to diffuse to the center.

```cpp
for (auto&amp; [name, info]: rods){
    int x = atoi(name.substr(0,2).c_str());
    int y = atoi(name.substr(3,2).c_str());
    double my_amount = produced_n[name];

    for (const Direction&amp; dir : directions){
      int nx = x + dir.x; std::string nxs = "";
      int ny = y + dir.y; std::string nys = "";
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
```

Rods on the edge that has their neutrons diffusing into "empty space" will diffuse into the metal walls of the reactor, and does not affect fission at all, so we just remove it from the population all together.

```cpp
else{
        double transport_rate = k_factors[name];

        double lost_flux = my_amount * 0.94;
        double transported = (my_amount - lost_flux)/10 * transport_rate;

        diffusion_deltas[name] -= transported;
}
```

Then, to maintain proper temperature, heat must also be removed due to the boiling water converting into steam and cooling down, and heat that is removed from recirculating the water around the core. This is done using this bit of code. This section also simulate decay heat using the Wigner-Way formula formula for decay heat.

```cpp
  double time_since_sd = 0;
  double power_before_sd = 0;

  if (avg_power &lt; 0.02){
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

  double boilTemp = 110; //so much for the core being slightly pressurized smh
  if (model.reactorWaterTemperature &gt; boilTemp){
    double excess = model.reactorWaterTemperature - boilTemp;

    double steam_generated = excess * 1;
    model.reactorWaterTemperature -= steam_generated * 0.01;
  }

  model.reactorWaterTemperature += new_temp - old_temp;
```

Finally, the simulation can be ran by simply just running reactor_physics() with the delta time and the rods variable

```cpp
void setup() {
  Serial.begin(9600);
  for (const char* rod : rodPositions){
    Serial.println(rod);
    rods[std::string(rod)] = {68200531029,48};
  };
}

void loop() {
  float dt = 5.00/60.00;
  //temp rod controls
  int joystickValue = analogRead(34);
  if (joystickValue &gt;= 0 &amp;&amp; joystickValue &lt;= 2000){
    for (auto&amp; [name,info] : rods){
      info.p += 1;
      if (info.p &gt; 48){
        info.p = 48;
      }
    }
  }
  if (joystickValue &gt; 3500) {
    for (auto&amp; [name,info] : rods){
      info.p -= 1;
      if (info.p &lt; 0){
        info.p = 0;
      }
    }
  }
  if (digitalRead(2) == HIGH){
    for (auto&amp; [name,info] : rods){
      info.p = 0;
    }
  }
  reactor_physics(dt,rods);
  unsigned long long neutrons = 0;
  for (auto&amp; [name,info] : rods){
    neutrons += info.n;
  }
  double avg_flux = neutrons/rods.size();
  double APRM = (avg_flux/FullPowerFlux)*100;
  Serial.println(String(neutrons) + " C/s | " + String(APRM) + " %NOM | " + String(model.reactorWaterTemperature) + " Deg C | " + rods["15-03"].p + " Notches");
  sleep(dt);
}
```

This section also handles the temporary SCRAM (shutdown) switch and a joystick to allow for manual control rod movement, hence its length. Realistically, you would only need:

```cpp
void loop(){
  float dt = 5.00/60.00;
  reactor_physics(dt,rods);
  sleep(dt);
}
```

to simulate the reactor itself. This makes controls easier to hook up to the simulator as you would only have to change the rod positions which coudl easily be done by just setting any rod to an integer, for example:

```cpp
rods["15-03"].p = 4;
```

![Initiating a manual reactor SCRAM](./readme%20images/IMG_4347.gif)

A better version of this video can be viewed on [youtube](https://www.youtube.com/watch?v=QIrbZ2ZCz5E)