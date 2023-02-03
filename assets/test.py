import subprocess
import argparse

parser = argparse.ArgumentParser(description ='vina commands')
parser.add_argument('-r','--receptor',type=str, help='Receptor for processing')
parser.add_argument('-c','--center',type=str,required=True,help='center position')
parser.add_argument('-s','--size',type=str,required=True,help='box size')
parser.add_argument('-m','--module',type=str,required=True,help='module for docking')
parser.add_argument('-d','--docking',type=str,required=True,help='basic or flexible docking')

args = parser.parse_args()
receptor = args.receptor
module = args.module
docking = args.docking
center = (args.center).split(",")
center_x = float(center[0])
center_y = float(center[1])
center_z = float(center[2])
box = (args.size).split(",")
size_x = float(box[0])
size_y = float(box[1])
size_z = float(box[2])
#subprocess.run(["cat config.txt"],shell=True)

def readInput():
    '''    
       if True:
       print("This is an incorrect input file")
       quit()
    '''
    print(f'{receptor} {module} {center_x} {size_z}')
    f = open(f'{receptor}',"r")
    s = f.read()
    print(s)
#    lines = s.splitlines()
#    print(lines)

#    receptor = lines[0]
#    xcord = lines[1]
#    ycord = lines[2]
#    zcord = lines[3]

#    print(f'{receptor} {xcord} {ycord} {zcord}')

def main():
    readInput()

main()
