runDir = pwd()
cd("..")
homeDir = pwd()
if isfile(".userSettings.conf")
    mv(".userSettings.conf", ".userSettings.conf.bkp", remove_destination=true)
end
cp("test/.userSettings.conf", ".userSettings.conf", remove_destination=true)

testDir = joinpath(homeDir, "test")
addDataDir = joinpath(homeDir, "additionalData")
try
  run(`rm -r test/lib`)
catch
end
try
  run(`rm -r test/input`)
catch
end
try
  run(`rm -r test/output`)
catch
end
run(`julia Config.jl --tmpDir $testDir`)

kernelFile = joinpath(addDataDir, "testKernels.tm")
run(`julia Config.jl --kernelFile $kernelFile`)

meshFile = joinpath(addDataDir, "SHAP5.ply")
run(`julia Config.jl --meshFile $meshFile`)

meshFileShadow = joinpath(addDataDir, "SHAP5_shadow.ply")
run(`julia Config.jl --meshFileShadow $meshFileShadow`)

dataFile = joinpath(addDataDir, "gas.dat")
run(`julia Config.jl --dataFile $dataFile`)

# run test with gas case and no shadow checking
run(`julia Config.jl --doCheckShadow no`)
cd("src")
run(`julia main.jl 2015-02-11T00:00:00 TEST`)


# run test with gas case and with shadow checking
cd("..")
run(`julia Config.jl --doCheckShadow yes`)
cd("src")
run(`julia main.jl 2015-02-11T00:00:00 TEST`)


# run test with dust case and with shadow
dataFile = joinpath(addDataDir, "dust.dat")
cd("..")
run(`julia Config.jl --dataFile $dataFile`)
cd("src")
run(`julia main.jl 2015-02-11T00:00:00 TEST`)

# run test with dust case and no shadow
cd("..")
run(`julia Config.jl --doCheckShadow yes`)
cd("src")
run(`julia main.jl 2015-02-11T00:00:00 TEST`)
