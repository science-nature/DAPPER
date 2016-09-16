

import sys, os
menu_actions  = {}  
class bcolors:
#to color the instructions in the console
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'

# Main menu
def main_menu():
    os.system('clear')
    

    choice = raw_input(" >>  ")
    exec_menu(choice)
 
    return
ans=True
print(bcolors.OKBLUE +"""
Welcome,\n
Please choose the menu you want to start:
Chose which benchmark you want to launch
\n
1. Strongly coupled - fully observed - changing N
2. Strongly coupled - atmosphere observed - changing N
3. Strongly coupled - ocean observed - changing N
\n
4. Strongly coupled - fully observed - changing dtObs^full
5. Strongly coupled - atmosphere observed - changing dtObs^full
6. Strongly coupled - ocean observed - changing dtObs^full
\n
7. Strongly coupled - dtObs^oc = 2 weeks - changing dtobs^atm
8. Strongly coupled - dtObs^atm = 12 hours - changing dtobs^oc
\n
9. Weakly coupled - dtObs^oc = 2 weeks - changing dtobs^atm
10. Weakly coupled - dtObs^atm = 12 hours - changing dtobs^oc
\n
11. Quit
"""+ bcolors.ENDC)
ans=input(bcolors.OKGREEN +"What would you like to do? "+ bcolors.ENDC)
if ans=="1": 
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/obsFULL"+ bcolors.ENDC)
  from benchmarks.benchmark2 import *
elif ans=="2":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/obsATM"+ bcolors.ENDC)
  from benchmarks.benchmark2atm import *
elif ans=="3":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/obsOC"+ bcolors.ENDC)
  from benchmarks.benchmark2oc import *
elif ans=="4":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/dtObs"+ bcolors.ENDC)
  from benchmarks.benchmark3 import *
elif ans=="5":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/obsATMdtObs"+ bcolors.ENDC)
  from benchmarks.benchmark3atm import *
elif ans=="6":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/obsOCdtObs"+ bcolors.ENDC)
  from benchmarks.benchmark3oc import *
elif ans=="7":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/dtobsATM"+ bcolors.ENDC)
  from benchmarks.benchmark3dtatm import *
elif ans=="8":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/strongly/dtobsOC"+ bcolors.ENDC)
  from benchmarks.benchmark3dtoc import *
elif ans=="9":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/weakly/dtobsATM"+ bcolors.ENDC)
  from benchmarks.benchmark3watm import *
elif ans=="10":
  print(bcolors.OKBLUE +"You will find the data in the file ./data/weakly/dtobsOC"+ bcolors.ENDC)
  from benchmarks.benchmark3woc import *
elif ans=="11":
  exit()
