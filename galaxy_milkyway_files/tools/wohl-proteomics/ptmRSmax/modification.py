class Modification:
    def __init__(self):
        self.AAs=[]
        self.mass=0.0
        self.NL={}
        self.DI={}
        self.factor=1
        self.name=""
        self.ID=-1

    def setAAs(self, newAAlist):
        self.AAs=newAAlist
        return True

    def getAAs(self):
        return self.AAs

    def setMass(self, newMass):
        self.mass=newMass
        return True

    def getMass(self):
        return self.mass

    def setNL(self, newNLdict):
        self.NL=newNLdict
        return True

    def getNL(self):
        return self.NL

    def setDI(self, newDIdict):
        self.DI=newDIdict
        return True

    def getDI(self):
        return self.DI

    def setName(self, newName):
        self.name=newName
        return True

    def getName(self):
        return self.name

    def setID(self, newID):
        self.ID=str(newID)
        return True

    def getID(self):
        return self.ID

    def create(self,nameString,modString,NLstring):#,DIstring):
        self.name=nameString
        
        for each in modString.split()[0]: # This will be each AA...
            self.AAs.append(each)
        self.mass=modString.split()[1]

        for each in NLstring.split()[0]:
            if each != "-":
                self.NL[each]=float(NLstring.split()[1])
        
        #for each in DIstring.split()[0]:
        #    if each != "-":
        #        self.DI[each]=float(DIstring.split()[1])

    def __str__(self):
        return self.ID+"_"+self.name+"_"+str(self.mass)
