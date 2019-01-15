import inLib
import time
import numpy as np
import subprocess
from scipy.ndimage import rotate

class Control(inLib.Device):
    def __init__(self, settings):

        self.pixels = settings['pixels']
        self.segments = settings['segments']
        self.executable = settings['executable']

        self.tempfilename = 'mirrorSegs.txt'
        self.multiplier = 1.0

        self.preMultiplier = 80

        self.zernike = None


        self.group = []

        self.padding = False

        self.proc = None
        self.numZernsToVary = 100

    
    def setup_dummy_mirror(self):
        pass

    def highlight_dummy_mirror_segs(self, segs):
        self.group = segs
        for i in range(len(segs)):
            self.dummy_mirror.addOffset(segs[i], 1)
        return self.dummy_mirror.pattern

    def clearDummyMirror(self):
        self.dummy_mirror.clear()

    def returnDummyMirror(self):
        return self.dummy_mirror.pattern, self.dummy_mirror.segOffsets
        

    def loadPattern(self, pattern_filename, mult=1.0):
        self.mirror.pattern = np.load(str(pattern_filename)) * mult
        return self.mirror.pattern

    def loadSegments(self, filename):
        msegs = np.loadtxt(filename)
        #print "msegs: ", msegs
        self.mirror.inputMirrorSegs(msegs)
        
    def patternRot90(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.rot90(self.mirror.pattern)
        return self.mirror.pattern

    def patternFlipLR(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.fliplr(self.mirror.pattern)
        return self.mirror.pattern

    def patternFlipUD(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.flipud(self.mirror.pattern)
        return self.mirror.pattern

    def patternRotate(self, deg):
        if self.mirror.pattern is not None:
            self.mirror.pattern = rotate(self.mirror.pattern, deg)
        return self.mirror.pattern

    def pokeGroup(self, group, offset, quiet=False):
        for s in group:
            self.pokeSegment(s, offset)
            
    def pokeSegment(self, seg, value, pokeAll = False):
        '''
        Add value to a given segment or to all segments
        '''
        if pokeAll:
            self.mirror.addOffset(-1,value)
        else:
            self.mirror.addOffset(seg,value)

    def reconfigGeo(self, cx, cy, npixels):
        self.mirror.initGeo(npixels)
        half_px = npixels/2
        if self.mirror.pattern is not None:
            if self.mirror.pattern.shape[0] > npixels:
                pattern = self.mirror.pattern[cx-half_px:cx+half_px,
                                              cy-half_px:cy+half_px]
                if (pattern.shape[0] == npixels) and (pattern.shape[1]==npixels):
                    self.mirror.pattern = pattern.copy()
                else:
                    print("Problem cropping pattern for DM.")
            if self.mirror.pattern.shape[0] < npixels:
                pattern = np.zeros((npixels,npixels))
                origShape = self.mirror.pattern.shape[0]
                pattern[cx-origShape/2:cx+origShape/2, cy-origShape/2:cy+origShape/2] = self.mirror.pattern
                self.mirror.pattern = pattern.copy()
        self._geometry = self.mirror.geometry
        return self.returnPattern()

    def padZeros(self, border, always=False):
        '''
        Pad pattern with zeros
        '''
        currentPixels = self.mirror.pattern.shape[0]
        newPattern = np.zeros((currentPixels+2*border, currentPixels+2*border))
        newPattern[border:-1*border,border:-1*border] = self.mirror.pattern
        self.mirror.initGeo(newPattern.shape[0])
        self.mirror.setPattern(newPattern)
        return self.returnPattern()

    def findSegments(self):
        self.mirror.findSegOffsets()

    def setMultiplier(self,mult):
        self.multiplier = mult

    def setPreMultiplier(self,mult):
        self.preMultiplier = mult

    def getSegments(self):
        return self.mirror.returnSegs()
        #return self.mirror.segOffsets

    def returnSegments(self):
        return self.mirror.returnSegs()
        #return self.mirror.segOffsets

    def returnPattern(self):
        return self.mirror.pattern

    def clear(self):
        self.mirror.clear()

    def setZernMode(self, mode):
        self.zernMode = mode
    
    def calcZernike(self, mode, amp, radius=None, useMask=True):
        if radius is None:
            radius = self.mirror.nPixels/2
        modes = np.zeros((mode))
        modes[mode-1]=amp
        self.zernike = libtim.zern.calc_zernike(modes, radius, mask=useMask,
                                                zern_data = {})
        return self.zernike

    def addZernike(self, zernike_pattern=None):
        if zernike_pattern is not None:
            zern = zernike_pattern
        else:
            if self.zernike is None:
                return 0
            else:
                zern = self.zernike
        if self.zernike is not None:
            p=self.mirror.addToPattern(zern * self.preMultiplier)
        return self.returnPattern()

    def getNumberOfZernToVary(self):
        return self.numZernsToVary

    def varyMultiplierCurrent(self, minMult, maxMult, num, wTime, externallyCalled = False):
        print(("The number of multipliers:", num))
        print(("mininum:", minMult))
        print(("maximum:",maxMult))
        self.numZernsToVary = num
        mults = np.linspace(minMult,maxMult,num)
        filenms = []
        baseline = self.returnPattern()
        for i in range(num):
            self.clear()
            newPattern = baseline * mults[i]
            self.mirror.setPattern(newPattern)
            self.findSegments()
            filenms.append("segFile_mult_%.2i.txt" % i)
            self.mirror.outputSegs(filenms[i])
        files_file = "allMultFiles_Max%.2i.txt" % maxMult
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print("Finished creating files for varying multiplier for current pattern...")

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr]

        if self.proc is not None:
            print("Polling proc: ", self.proc.poll())
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print("Going to next...")
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print("proc stdout: ", output)
            return 1
        

    def varyZernAmp(self, mode, maxAmp, minAmp, num, wTime, radius=None, useMask=True, clearfirst=True,
                    externallyCalled = False):
        '''
        Calls external C++ program that applies *num* different patterns of
        Zernike mode *mode* to the mirror.

        :Parameters:
            *mode*: int
                Zernike mode
            *maxAmp*: float
            *minAmp*: float
            *num*: int
            *wTime*: float
                Time to wait in milliseconds before new pattern applied to mirror
            *useMask*: boolean
                Optional
            *clearFirst*: boolean
                Optional
        '''
        #mode = self.zernMode
        self.numZernsToVary = num
        amps = np.linspace(minAmp,maxAmp,num)
        filenms = []
        baseline = self.returnPattern()
        for i in range(num):
            if clearfirst:
                self.clear() #clears mirror pattern and segments
            else:
                self.clear()
                self.mirror.setPattern(baseline)
            zern = self.calcZernike(mode, amps[i], radius=radius, useMask=useMask)
            self.addZernike(zernike_pattern=zern)
            self.findSegments()
            filenms.append("segFile_mode%.3i_amp%.2i.txt" % (mode,i))
            self.mirror.outputSegs(filenms[i])
        files_file = "allSegFiles_%.3i.txt" % mode
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print("Finished creating files for varying Zernike...")

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr] 

        #subprocess.call([self.executable, files_file, str(self.multiplier),
        #                 numStr, wTimeStr], shell=True)

        if self.proc is not None:
            print("Polling proc: ", self.proc.poll())
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print("Going to next...")
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print("proc stdout: ", output)
            return 1

    def varyZernRadii(self, mode, amp, maxR, minR, num, wTime, radius=None, useMask=True, clearfirst=True,
                      externallyCalled = False):
        '''
        Calls external C++ program that applies *num* different patterns of
        Zernike mode *mode* to the mirror.

        :Parameters:
            *mode*: int
                Zernike mode
            *amp*: float
                Zernike amplitude
            *maxR*: float
            *minR*: float
            *num*: int
            *wTime*: float
                Time to wait in milliseconds before new pattern applied to mirror
            *useMask*: boolean
                Optional
            *clearFirst*: boolean
                Optional
        '''
        #mode = self.zernMode
        self.numZernsToVary = num
        rads = np.linspace(minR,maxR,num,dtype=np.uint16)
        filenms = []
        baseline = self.returnPattern()
        for i in range(num):
            if clearfirst:
                self.clear() #clears mirror pattern and segments
            else:
                self.clear()
                self.mirror.setPattern(baseline)
            zern = self.calcZernike(mode, amp, radius=rads[i], useMask=useMask)
            self.addZernike(zernike_pattern=zern)
            self.findSegments()
            filenms.append("segFile_mode%.3i_rad%.2i.txt" % (mode,i))
            self.mirror.outputSegs(filenms[i])
        files_file = "allSegFiles_%.3i.txt" % mode
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print("Finished creating files for varying Zernike radii...")

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr] 

        #subprocess.call([self.executable, files_file, str(self.multiplier),
        #                 numStr, wTimeStr], shell=True)

        if self.proc is not None:
            print("Polling proc: ", self.proc.poll())
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print("Going to next...")
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print("proc stdout: ", output)
            return 1
        

    def advancePatternWithPipe(self):
        if self.proc is not None:
            self.proc.stdin.write("\n")
            
            
            
    def advancePipe(self):
        print("The process:", self.proc)
        if self.proc is not None:
            print(self.proc.communicate())
            self.proc = None

   
    
    def applyToMirror(self, wtime=-1):
        #First save mirror
        self.mirror.outputSegs(self.tempfilename)

        #Wait to make sure file exists
        time.sleep(0.5)

        wTimeStr = str(wtime)

        self.proc = subprocess.Popen([self.executable, self.tempfilename, str(self.multiplier),"1", wTimeStr] , stdin = subprocess.PIPE, stdout = subprocess.PIPE)
        for ir in range(11):
            line = self.proc.stdout.readline()
            dec_line = line.rstrip().decode()
            print(dec_line)
            if line == '':
                print('Finished reading.')
                break
        print("The pattern is added to the mirror.")
        #Run executable
        #subprocess.call([self.executable, self.tempfilename, str(self.multiplier),"1", wTimeStr], shell=True)
        #subprocess.call([self.executable, self.tempfilename, str(self.multiplier),
         #                "1", "30000"], shell=True)
    