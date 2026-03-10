local function wait(seconds)
    local start = os.clock()
    repeat
    until os.clock() > start + seconds
end

local function clamp(val, min, max)
    return math.max(min, math.min(max, val))
end

local rods = { --where n is neutrons, p is rod position (0 is fully inserted, 48 is fully withdrawn)
    ["15-03"] = {n = 0, p = 0}, ["21-03"] = {n = 0, p = 0}, ["27-03"] = {n = 0, p = 0},
    ["09-09"] = {n = 0, p = 0}, ["15-09"] = {n = 0, p = 0}, ["21-09"] = {n = 0, p = 0}, ["27-09"] = {n = 0, p = 0}, ["33-09"] = {n = 0, p = 0},
    ["03-15"] = {n = 0, p = 0}, ["09-15"] = {n = 0, p = 0}, ["15-15"] = {n = 0, p = 0}, ["21-15"] = {n = 0, p = 0}, ["27-15"] = {n = 0, p = 0}, ["33-15"] = {n = 0, p = 0}, ["39-15"] = {n = 0, p = 0},
    ["03-21"] = {n = 0, p = 0}, ["09-21"] = {n = 0, p = 0}, ["15-21"] = {n = 0, p = 0}, ["21-21"] = {n = 0, p = 0}, ["27-21"] = {n = 0, p = 0}, ["33-21"] = {n = 0, p = 0}, ["39-21"] = {n = 0, p = 0},
    ["03-27"] = {n = 0, p = 0}, ["09-27"] = {n = 0, p = 0}, ["15-27"] = {n = 0, p = 0}, ["21-27"] = {n = 0, p = 0}, ["27-27"] = {n = 0, p = 0}, ["33-27"] = {n = 0, p = 0}, ["39-27"] = {n = 0, p = 0},
    ["09-33"] = {n = 0, p = 0}, ["15-33"] = {n = 0, p = 0}, ["21-33"] = {n = 0, p = 0}, ["27-33"] = {n = 0, p = 0}, ["33-33"] = {n = 0, p = 0},
    ["15-39"] = {n = 0, p = 0}, ["21-39"] = {n = 0, p = 0}, ["27-39"] = {n = 0, p = 0}
}


local reactor_inventory = {
    waterMass = 1000
}

local model = {
    reactor_water_temperature = 20,
}

local BARN = 1e-24
local N235 = 2.53e22
local N238 = 2.51e22

local xs = {
    U235_capture = 99,
    U235_fission = 583,
    U238_capture = 2,
    H2_scatter = 4,
    C12_capture = 0.002,
    B10_capture = 200
}

local function fuel_get(waterMass, controlDepth, neutronFlux, temperatureFuel, coreFlow)
    
    local Sigma_f_235 = N235 * xs.U235_fission * BARN
    local Sigma_c_235 = N235 * xs.U235_capture * BARN
    local Sigma_c_238 = N238 * xs.U238_capture * BARN
    
    local Sigma_a_fuel = Sigma_f_235 + Sigma_c_235 + Sigma_c_238

    local v = 2.43
    local eta = (v * Sigma_f_235) / Sigma_a_fuel

    local U = Sigma_a_fuel
    
    local M = xs.C12_capture * BARN * N235 
    
    local insertion_amount = (1 - controlDepth)
    local CR = (xs.B10_capture * BARN * N235) * (insertion_amount * 2.5)

    local f = U / (U + M + CR)

    local voids = neutronFlux / 2.5e12
    coreFlow = clamp(math.abs((coreFlow/100)-1),0,1)
    voids = clamp(coreFlow*voids,0,0.7)

    f = f*(1-voids)

    local p = 0.9
    local eps = 1.03
    local Pfnl = 0.97
    local Ptnl = 0.97

    local kEff = eta * f * p * eps * Pfnl * Ptnl

    local kStep = kEff^0.03

    return {
        kStep = kStep,
        kEff = kEff
    }
end

local function reactor_physics(dt, rods)
    local CoreFlow = 0
    local waterMass = reactor_inventory.waterMass
    
    local next_rods_n = {}

    for name, info in pairs(rods) do
        local NeutronFlux = math.max(info.n, 100)
        local mykEffArgs = fuel_get(waterMass, info.p/48, NeutronFlux, 60, CoreFlow)
        local mykStep = mykEffArgs.kStep
        
        next_rods_n[name] = math.max(info.n * mykStep, 10)
    end

    local directions = {
        {x=6,y=0},
        {x=-6,y=0},
        {x=0,y=6},
        {x=0,y=-6}
    }

    local diffusion_changes = {}

    for name, info in pairs(rods) do
        local x, y = name:match("(%d+)%-(%d+)")
        x = tonumber(x)
        y = tonumber(y)

        local current_n = rods[name].n

        for _, dir in pairs(directions) do
            local nx = x + dir.x
            local ny = y + dir.y
            local neighbor_name = string.format("%02d-%02d", nx, ny)

            if rods[neighbor_name] then
                local neighbor_n = rods[neighbor_name].n
                
                local mykEffArgs = fuel_get(waterMass, info.p/48, math.max(current_n, 100), 60, CoreFlow)
                local transport_rate = mykEffArgs.kStep 
                
                local transported = (current_n - neighbor_n) / 10 * transport_rate

                diffusion_changes[name] = (diffusion_changes[name] or 0) - transported
                diffusion_changes[neighbor_name] = (diffusion_changes[neighbor_name] or 0) + transported
            end
        end
    end

    local total_energy = 0
    
    for name, n_val in pairs(next_rods_n) do
        local diff = (diffusion_changes[name] or 0) / 2
        local final_n = math.max(n_val + diff, 10)
        
        rods[name].n = final_n
        
        total_energy = total_energy + (final_n / 2.5e12)
    end

    local energy = total_energy * 3486
    local calories = energy * 1000000
    local HeatC = calories / 1000
    local TempChange = (HeatC / math.max(waterMass, 1)) * dt
    model.reactor_water_temperature = model.reactor_water_temperature + TempChange
end

while true do
    local dt = 5/60
    reactor_physics(dt, rods)
    
    local total_n = 0
    for k,v in pairs(rods) do total_n = total_n + v.n end
    print("Total Neutrons: " .. total_n)
    print("Moderator Temperature: " .. model.reactor_water_temperature)
    
    wait(dt)
end