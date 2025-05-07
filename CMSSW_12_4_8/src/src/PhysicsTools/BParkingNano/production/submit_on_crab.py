from CRABClient.UserUtilities import config, ClientException
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'TTToHadronic{:s}'.format(production_tag)

config.section_('Data')
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/gcelotto/btv_ntuples/{:s}'.format(config.General.workArea)
config.Data.inputDBS = 'global'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/run_nano_cfg.py'
config.JobType.maxJobRuntimeMin = 2750
config.JobType.allowUndistributedCMSSW = True

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_CH_PSI'

if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from http.client import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config)
      except HTTPException as hte:
          print("Failed submitting task:",hte.headers)
      except ClientException as cle:
          print("Failed submitting task:",cle)

  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  parser.add_argument('-r', '--lhcRun', type=int, default=2, help = 'Run 2 or 3 (default)')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f,Loader=yaml.FullLoader) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}

    # loop over samples
    for sample, info in doc['samples'].items():
      # Input DBS
      input_dbs = info['dbs'] if 'dbs' in info else None
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample.replace('%d',str(part)) if part is not None else sample

        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print('submitting', name)

        isMC = info['isMC']

        config.Data.inputDBS = input_dbs if input_dbs is not None else 'global'

        config.Data.inputDataset = info['dataset'].replace('%d',str(part)) \
                                   if part is not None else \
                                      info['dataset']

        config.General.requestName = name
        common_branch = 'mc' if isMC else 'data'
        config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        if not isMC:
            config.Data.lumiMask = info.get(
                'lumimask',
                common[common_branch].get('lumimask', None)
            )
        else:
            config.Data.lumiMask = ''

        config.Data.totalUnits = 1000
        config.Data.unitsPerJob = info.get(
            'splitting',
            common[common_branch].get('splitting', None)
        )
        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )

        config.JobType.pyCfgParams = [
            'isMC={:.0f}'.format(int(isMC)),
            'reportEvery=1000',
            'tag={:s}'.format(production_tag),
            'globalTag={:s}'.format(globaltag),
            'lhcRun={:.0f}'.format(args.lhcRun),
        ]

        config.JobType.outputFiles = ['_'.join(['TTToHadronic', 'Run3' if args.lhcRun==3 else 'Run2', 'mc' if isMC else 'data', production_tag])+'.root']

        print()
        print(config)
        submit(config)
