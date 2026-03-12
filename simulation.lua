local function wait(seconds)
    local start = os.clock()
    repeat
    until os.clock() > start + seconds
end

local function clamp(val, min, max)
    return math.max(min, math.min(max, val))
end

function string:split(sep)
    local sep, fields = sep or " ", {}
    local pattern = string.format("([^%s]+)", sep)
    self:gsub(pattern, function(c) fields[#fields+1] = c end)
    return fields
end

local rods = { 
    
    ["15-03"] = {n = 0, p = 0}, ["21-03"] = {n = 0, p = 0}, ["27-03"] = {n = 0, p = 0},
    ["09-09"] = {n = 0, p = 0}, ["15-09"] = {n = 0, p = 0}, ["21-09"] = {n = 0, p = 0}, ["27-09"] = {n = 0, p = 0}, ["33-09"] = {n = 0, p = 0},
    ["03-15"] = {n = 0, p = 0}, ["09-15"] = {n = 0, p = 0}, ["15-15"] = {n = 0, p = 0}, ["21-15"] = {n = 0, p = 0}, ["27-15"] = {n = 0, p = 0}, ["33-15"] = {n = 0, p = 0}, ["39-15"] = {n = 0, p = 0},
    ["03-21"] = {n = 0, p = 0}, ["09-21"] = {n = 0, p = 0}, ["15-21"] = {n = 0, p = 0}, ["21-21"] = {n = 0, p = 0}, ["27-21"] = {n = 0, p = 0}, ["33-21"] = {n = 0, p = 0}, ["39-21"] = {n = 0, p = 0},
    ["03-27"] = {n = 0, p = 0}, ["09-27"] = {n = 0, p = 0}, ["15-27"] = {n = 0, p = 0}, ["21-27"] = {n = 0, p = 0}, ["27-27"] = {n = 0, p = 0}, ["33-27"] = {n = 0, p = 0}, ["39-27"] = {n = 0, p = 0},
    ["09-33"] = {n = 0, p = 0}, ["15-33"] = {n = 0, p = 0}, ["21-33"] = {n = 0, p = 0}, ["27-33"] = {n = 0, p = 0}, ["33-33"] = {n = 0, p = 0},
    ["15-39"] = {n = 0, p = 0}, ["21-39"] = {n = 0, p = 0}, ["27-39"] = {n = 0, p = 0}
}

local reactor_inventory = {
    waterMass = 75000 
}

local model = {
    reactor_water_temperature = 20, 
    APRM_AVG = 0, 
    IRM = {}, 
    irm_range_switch = 1 
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
    B10_capture = 3840
}

local FULL_POWER_FLUX = 2.5e12
local IRM_SCALES = {}

local ir_top_scale_percent = 40    
local ir_bottom_scale_percent = 0.00004 

for i = 1, 12 do
    
    local pct = ir_bottom_scale_percent * (ir_top_scale_percent / ir_bottom_scale_percent) ^ ((i-1)/11)
    IRM_SCALES[i] = FULL_POWER_FLUX * (pct / 100)
end

local function fuel_get(waterMass, controlDepth, neutronFlux, temperatureFuel, coreFlow)
    local Sigma_f_235 = N235 * xs.U235_fission * BARN
    local Sigma_c_235 = N235 * xs.U235_capture * BARN
    local Sigma_c_238 = N238 * xs.U238_capture * BARN
    
    local Sigma_a_fuel = Sigma_f_235 + Sigma_c_235 + Sigma_c_238

    local v = 2.43
    local eta = (v * Sigma_f_235) / Sigma_a_fuel

    local U = Sigma_a_fuel
    local M = xs.C12_capture * BARN * N235 
    
    local insertion = 1 - controlDepth 
    local CR = (xs.B10_capture * BARN * N235) * insertion 

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

    return { kStep = kStep, kEff = kEff }
end

local function reactor_physics(dt, rods)
    local CoreFlow = 0
    local waterMass = reactor_inventory.waterMass
    
    local produced_n = {}
    local k_factors = {} 

    for name, info in pairs(rods) do
        local NeutronFlux = math.max(info.n, 100)
        local mykEffArgs = fuel_get(waterMass, info.p/48, NeutronFlux, 60, CoreFlow)
        
        k_factors[name] = mykEffArgs.kStep
        produced_n[name] = math.max(info.n * mykEffArgs.kStep, 10)
    end
    
    local diffusion_deltas = {}
    local directions = {{x=4,y=0}, {x=-4,y=0}, {x=0,y=4}, {x=0,y=-4}}

    for name, _ in pairs(rods) do
        local x, y = name:match("(%d+)%-(%d+)")
        x = tonumber(x)
        y = tonumber(y)
        local my_amount = produced_n[name]
        
        for _, dir in pairs(directions) do
            local nx = x + dir.x
            local ny = y + dir.y
            local neighbor_name = string.format("%02d-%02d", nx, ny)

            if produced_n[neighbor_name] then
                local neighbor_amount = produced_n[neighbor_name]
                local transport_rate = k_factors[name] 
                local transported = (my_amount - neighbor_amount) / 10 * transport_rate
                
                diffusion_deltas[name] = (diffusion_deltas[name] or 0) - transported
                diffusion_deltas[neighbor_name] = (diffusion_deltas[neighbor_name] or 0) + transported
            end
        end
    end

    
    local total_neutrons = 0
    local total_lprm_sum = 0
    local rod_count = 0
    
    for name, val in pairs(produced_n) do
        local diff = (diffusion_deltas[name] or 0) / 2
        local final = math.max(val + diff, 10)
        
        rods[name].n = final
        total_neutrons = total_neutrons + final
        
        local srm_val = math.min(final * 10, 1e6)
        rods[name].srm = srm_val
        
        local current_range = model.irm_range_switch or 1
        
        local scale_max_flux = IRM_SCALES[current_range] or IRM_SCALES[12]
        
        local irm_deflection = (final / scale_max_flux) * 125
        rods[name].irm = clamp(irm_deflection, 0, 125)
        
        local lprm_val = (final / FULL_POWER_FLUX) * 100
        rods[name].lprm = lprm_val
        rods[name].aprm = lprm_val
        
        total_lprm_sum = total_lprm_sum + lprm_val
        rod_count = rod_count + 1
    end

    if rod_count > 0 then
        model.APRM_AVG = total_lprm_sum / rod_count
    end
    
    local power_fraction = model.APRM_AVG / 100
    local current_mwt = power_fraction * 760 --185 rods = 3800MWt, so scaled down to 37 rods is 760MWt (adjustable)
    local kcal_per_sec = current_mwt * 239
    
    local temp_change = (kcal_per_sec / waterMass) * dt --todo: fix thermal hydraulics
    
    model.reactor_water_temperature = model.reactor_water_temperature + temp_change
end

local function process_command()
    io.write("\nCOMMAND (scram, set <id> <0-48>, range <1-12>, run <sec>, status) > ")
    local input = io.read()
    
    if not input or input == "" then return end
    
    local args = input:split(" ")
    local cmd = args[1]:lower()

    if cmd == "scram" then
        for k, v in pairs(rods) do v.p = 0 end

    elseif cmd == "set" then
        local target = args[2]
        local val = tonumber(args[3])
        if target == "all" then
            for k, v in pairs(rods) do v.p = val end
            print("All rods set to " .. val)
        elseif rods[target] then
            rods[target].p = val
            print(target .. " set to " .. val)
        else
            print("Invalid target.")
        end

    elseif cmd == "range" then
        local r = tonumber(args[2])
        if r and r >= 1 and r <= 12 then
            model.irm_range_switch = r
            
            local range_flux = IRM_SCALES[r]
            local range_pct = (range_flux / FULL_POWER_FLUX) * 100
            print(string.format("IRM Range %d selected (Max Scale: %.4f%% Power)", r, range_pct))
        else
            print("Invalid Range. Use 1-12.")
        end

    elseif cmd == "run" then
        local duration = tonumber(args[2]) or 1
        local elapsed = 0
        local dt = 5/60 
        
        print(string.format("Simulating for %.1f seconds...", duration))
        
        while elapsed < duration do
            reactor_physics(dt, rods)
            elapsed = elapsed + dt
            
            if math.fmod(elapsed, 0.5) < dt then
                print(string.format("T_POskut%.2f | N_SR:%.1f C/s | N_IR:%.1f %% | N_RE:%.0f %%", model.reactor_water_temperature, rods["21-21"].srm,rods["21-21"].lprm,model.APRM_AVG))
    
            end
            wait(dt)
        end

    elseif cmd == "status" then
        local count = 0
        for k, v in pairs(rods) do
            io.write(string.format("%s:%02d ", k, v.p))
            count = count + 1
            if count % 6 == 0 then io.write("\n") end
        end
        io.write("\n")

    elseif cmd == "exit" or cmd == "quit" then
        os.exit()
    end
end

while true do
    process_command()
end