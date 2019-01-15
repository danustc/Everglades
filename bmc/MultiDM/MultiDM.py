import inLib
import time
import numpy as np
import subprocess
from scipy.ndimage import rotate
from myWidget.DM_simulate import DM

class Control(inLib.Device):
    def __init__(self, settings):

        self.pixels = settings['pixels']
        self.segments = settings['segments']
        self.executable = settings['executable']

        self.tempfilename = 'mirrorSegs.txt'
        self.multiplier = 1.0

        self.preMultiplier = 80
        self.zernike = None
        self.mirror = DM(self.segments, self.pixels)
        self.mirror.findSeg()
        self.group = []

        self.proc = None



    def highlight_dummy_mirror_segs(self, segs):
        self.group = segs
        for i in range(len(segs)):
            self.dummy_mirror.addOffset(segs[i], 1)
        return self.dummy_mirror.pattern

    def clearDummyMirror(self):
        self.dummy_mirror.clear()

    def returnDummyMirror(self):
        return self.dummy_mirror.pattern, self.dummy_mirror.segOffsets
        
    def setPattern(self, newPattern):
        '''
        pass the new pattern to the Deformable mirror
        '''
        print("I received the new pattern!")


    def loadPattern(self, pattern_filename, mult=1.0):
        self.mirror.pattern = np.load(str(pattern_filename)) * mult
        return self.mirror.pattern

    def loadSegments(self, filename):
        msegs = np.loadtxt(filename)
        #print "msegs: ", msegs
        self.mirror.readSeg(msegs)
        
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
        self.mirror.clearPattern()

    def setZernMode(self, mode):
        self.zernMode = mode



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
        self.mirror.exportSegs(self.tempfilename)

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
    