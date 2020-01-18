import os
import sys
import shutil

if any(i in ['mercury_1', 'mercury_2', 'mercury_3', 'mercury_4'] for i in os.listdir('.')):
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n',
          'make.py:\n',
          'Four instances of mercury already exist. Exiting make.py.\n',
          '~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.exit()
elif 'mercury' not in os.listdir('.'):
    print(' ~~~~~~~~~~~~~~~~~~~~~~~~\n',
          'make.py:\n',
          'Original mercury/ directory not present. Exiting make.py.\n',
          '~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.exit()
else:
    # Compile files in main folder
    os.system('gfortran -o mercury/mercury.exe mercury/mercury6_2_mig.for')
    os.system('gfortran -o mercury/element.exe mercury/element6.for')
    os.system('gfortran -o mercury/close.exe mercury/close6.for')

    # Create 4 directories and copy over
    pnos = [1, 2, 3, 4]
    for pno in pnos:
        os.mkdir('mercury_{}'.format(pno))
        shutil.copyfile('mercury/README.txt', 'mercury_{}/README.txt'.format(pno))
        shutil.copyfile('mercury/big.in', 'mercury_{}/big.in'.format(pno))
        shutil.copyfile('mercury/close.in', 'mercury_{}/close.in'.format(pno))
        shutil.copyfile('mercury/close6.for', 'mercury_{}/close6.for'.format(pno))
        shutil.copyfile('mercury/element.in', 'mercury_{}/element.in'.format(pno))
        shutil.copyfile('mercury/element6.for', 'mercury_{}/element6.for'.format(pno))
        shutil.copyfile('mercury/files.in', 'mercury_{}/files.in'.format(pno))
        shutil.copyfile('mercury/mercury.inc', 'mercury_{}/mercury.inc'.format(pno))
        shutil.copyfile('mercury/mercury6.man', 'mercury_{}/mercury6.man'.format(pno))
        shutil.copyfile('mercury/mercury6_2_mig.for', 'mercury_{}/mercury6_2_mig.for'.format(pno))
        shutil.copyfile('mercury/message.in', 'mercury_{}/message.in'.format(pno))
        shutil.copyfile('mercury/param.in', 'mercury_{}/param.in'.format(pno))
        shutil.copyfile('mercury/small.in', 'mercury_{}/small.in'.format(pno))
        shutil.copyfile('mercury/swift.inc', 'mercury_{}/swift.inc'.format(pno))
        shutil.copyfile('mercury/mercury.exe', 'mercury_{0}/mercury_{0}.exe'.format(pno))
        shutil.copyfile('mercury/element.exe', 'mercury_{0}/element_{0}.exe'.format(pno))
        shutil.copyfile('mercury/close.exe', 'mercury_{0}/close_{0}.exe'.format(pno))

        # Make executable
        os.system('chmod +x mercury_{0}/mercury_{0}.exe'.format(pno))
        os.system('chmod +x mercury_{0}/element_{0}.exe'.format(pno))
        os.system('chmod +x mercury_{0}/close_{0}.exe'.format(pno))

