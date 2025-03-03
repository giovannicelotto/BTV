pyFile=TOP-RunIISummer20UL18MiniAODv2-00149_1_cfg.py
evContent=MINIAODSIM
custom=Configuration/DataProcessing/Utils.addMonitoring
datatier=MINIAOD
EVENTS=4200
NAME="ttbarFullyHadronic"
inFile="/store/mc/RunIISummer20UL18MiniAODv2/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/004EF875-ACBB-FE45-B86B-EAF83448CE62.root"
outFile=file:TOP-RunIISummer20UL18MiniAODv2-00149.root
cond=106X_upgrade2018_realistic_v16_L1v1

#step=RAW2DIGI,RECO,RECOSIM,PAT
step=PAT
procModifiers=run2_miniAOD_UL
geometry=DB:Extended
era=Run2_2018

# --inputCommands "keep *"
# --customise_commands "process.MINIAODSIMoutput.outputCommands.append('keep *_*_*_HLT'); process.MINIAODSIMoutput.outputCommands.append('keep *_*_*_SIM');"

cmsDriver.py  --python_filename $pyFile --eventcontent $evContent --customise $custom --datatier $datatier --fileout $outFile --conditions $cond --step $step --procModifiers $procModifiers --geometry $geometry --filein $inFile --era $era --runUnscheduled --no_exec --mc -n $EVENTS --dump_python || echo "Error" $? ;
