
import os
import random

try:
    os.mkdir('database')
except FileExistsError:
    pass

class Peptides():
    def __init__(self, initSeq, addResNum):        
        self.initSeq = initSeq
        self.addResNum = addResNum
        self.genNum = 1
        self.code_str = 'a r n d c q e g h i l k m f p s t w y v A R N D C Q E G H I L K M F P S T W Y V'
        self.code_list = self.code_str.split(" ")
        self.pep_dict = dict(zip(list(range(1, 41)), self.code_list))
        
    def genOneRandomPeptides(self):
        Seq = self.initSeq
        for j in range(1, self.addResNum+1):            
            Seq = Seq + self.pep_dict[random.randint(1,40)]
        return Seq
    
    def genMultipleRandomPeptides(self):
        f = open('./database/Peptides{}.txt'.format(str(len(self.initSeq) + self.addResNum))
            , mode='w'
            , encoding='utf-8'
        )
        for i in range(self.genNum):
            pep = self.genOneRandomPeptides()
            f.write(pep+'\n')
        f.close()
        
    def update_genNum(self, num):
        if num >= self.genNum:
            self.genNum = num
        else:
            print("Warning：You can't roll back an odometer!")
            
    def increment_genNum(self, num):
        self.genNum = self.genNum + num


my_pep = Peptides(initSeq = '', addResNum=6)
my_pep.genNum = 1000
my_pep.genMultipleRandomPeptides()

# my_pep.update_genNum(500)
# my_pep.genMultipleRandomPeptides()

# my_pep.increment_genNum(100)
# my_pep.genMultipleRandomPeptides()
# my_pep.increment_genNum(200)
# my_pep.genMultipleRandomPeptides()




