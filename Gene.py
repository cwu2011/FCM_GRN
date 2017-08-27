# -*- coding: utf-8 -*-
class Gene:
    def __init__(self):
        self.geneFamily_Str = ""
        self.geneMapped_Str = ""
        self.idGene_Str = ""
        self.expGene_List = []

    def setGene(self, geneFamily_Str, geneMapped_Str, idGene_Str, expGene_List):
        self.geneFamily_Str = geneFamily_Str
        self.geneMapped_Str = geneMapped_Str
        self.idGene_Str = idGene_Str
        self.expGene_List = list(expGene_List)

    def setGene0(self, idGene_Str, expGene_List):
        self.idGene_Str = idGene_Str
        self.expGene_List = list(expGene_List)
    def printGene(self):
        print str(self.geneFamily_Str)+' '+str(self.geneMapped_Str)+' '+str(self.idGene_Str)+' '+str(self.expGene_List)+'\n'