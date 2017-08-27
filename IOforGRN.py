# -*- coding: utf-8 -*-
import xlrd
import Gene

def readFirTimeSeriesUnigeneData(filename_Str):
    try:
        data = xlrd.open_workbook(filename_Str)
        table=data.sheets()[0] # sheets index starts from 0
        nrows=table.nrows #index starts from 0
        ncols=table.ncols #index starts from 0
        genes_List=[]
        for i in range(nrows):
            #print table.row_values(i)
            tr=table.row_values(i)
            genes_List.append(Gene.Gene())
            expList=[]
            for j in range (3 , ncols):
                expList.append(tr[j])
            genes_List[i].setGene(tr[0],tr[1],tr[2],expList)#tr[0-2]是基因的家族，名称，编号；tr[3-(ncols-1)(6 in this case)]是四个时间点下的表达量。如输入数据结构有变化，则须修改此处。
        # for i in range(nrows):
        #     genes_List[i].printGene()
        return genes_List
    except Exception, e:
        print str(e)

def readDREAM4Data(filanema_Str):

    return 0


# selcet lines need to annotate use the short cut: ctrl+/
# def main():
#     filename="E:\\data2read\GRN\\DataofFir\\"+"LOX_Family.xlsx"
#     genes_List=readFirTimeSeriesUnigeneData(filename)
#
# if __name__=="__main__":
#     main()


