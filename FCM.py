# -*- coding: utf-8 -*-
import numpy
from compiler.ast import flatten


######
# 计算目标函数 equ（5）
# function to call: fitness2(Args)
# parameters: c_2DArray(基因表达时间序列，标准化或sigmoid)；w_2DArray(权重矩阵，算法初始须生成一个)；lambda_sigmoid(sigmod的参数)

###############################

#####################################################

def get_CHats_2DArray(c_2DArray, w_2DArray):  # 2nd version of fitnessfunction
    c1_2DArray = c_2DArray[0:len(c_2DArray) - 1]  # 计算t+1时刻的估计，只需要前Tn-1个时间点的数据作为输入
    c1_Matrix = numpy.matrix(c1_2DArray)
    w_Matrix = numpy.matrix(w_2DArray)
    c_Hat_Matrix = c1_Matrix* w_Matrix
    c_Hat_2DArray = numpy.array(c_Hat_Matrix)
    # print c_Hat_2DArray
    return c_Hat_2DArray


def sigmoid_2DArray(lambda_sigmoid_Float, c_Har_2DArray):  # equ (4) g(·)转换函数，用于将活化度值限定在[0,1]内
    # ctt_array=numpy.array(ctt_matrix)
    # print ctt_array
    lam_x = -lambda_sigmoid_Float * c_Har_2DArray
    etemp = numpy.exp(lam_x)
    temp = 1.0 / (1.0 + etemp)  # type 2d array
    return temp


def fitness2(c_2DArray, w_2DArray, lambda_sigmoid_Float):
    # 是否需要对c_2DArray进行sigmoid计算？or 在作为参数之前就标准化或sigmoid化？
    # get C_HATS 不含最后一个时间点为Ct的计算结果
    c_Hat_2DArray = get_CHats_2DArray(c_2DArray, w_2DArray)
    # sigmoid
    c_Hat_2DArray_Sigmoid = sigmoid_2DArray(lambda_sigmoid_Float, c_Hat_2DArray)
    # 计算err
    data_error = 0
    m = (len(c_2DArray) - 1) * len(c_2DArray[0]) * 1.0  # (Nt-1)*Nn*Ns
    for i in range(1, len(c_2DArray)):  # t0 时刻的数据仅用于计算t1时刻的c_hat
        ct_array = c_2DArray[i]
        ct_hat_array = c_Hat_2DArray_Sigmoid[i - 1]
        delta = ct_array - ct_hat_array
        d_square = numpy.square(delta)
        d_sq_sum = numpy.sum(d_square)
        data_error += d_sq_sum
    print c_Hat_2DArray
    return data_error / m  # end of 2nd version of fitness function


############################
# def sigmoid(lambda_sigmoid, ctt_array):  # initial version of fitnessfunction # equ (4) g(·)转换函数，用于将活化度值限定在[0,1]内
#     # ctt_array=numpy.array(ctt_matrix)
#     # print ctt_array
#     lam_x = -lambda_sigmoid * ctt_array
#     etemp = numpy.exp(lam_x)
#     temp = 1.0 / (1.0 + etemp)  # type array
#     return temp
#
#
# def dFCM(w_array, ct_array, lambda_sigmoid):  # equ(3) caculate C(t+1)_hat using w and C(t)
#     # t_w_2Dlist = numpy.array(w_list).reshape(len(ct_list))
#     # w_matrix = numpy.mat(t_w_2Dlist)
#     w_matrix = numpy.mat(w_array)
#     c_matrix = numpy.mat(ct_array)
#     ctt_matrix = c_matrix * w_matrix
#     ctt_array = numpy.array(ctt_matrix)
#     # ctt_list = flatten(ctt_matrix.tolist())
#     return sigmoid(lambda_sigmoid, ctt_array)# C(t+1)_hat_array
#     #return ctt_array
#
#
# def getCHat(c_array, lambda_sigmoid, w_array):
#     nt = len(c_array)  # c_array 的行数即时间点的数量
#     c_hat_list = []
#     for i in range(1, nt):
#         c_hat_list.append(list(dFCM(w_array, c_array[i - 1], lambda_sigmoid)))
#     c_hat_array = numpy.array(c_hat_list)
#     return c_hat_array  #
#
#
# def objectiveFunction(c_array, c_hat_array):  # Ns=1 equ (5),c_list各行存储一个时间点下所有基因的表达量
#     data_error = 0
#     m = (len(c_array) - 1) * len(c_array[0]) * 1.0  # (Nt-1)*Nn*Ns
#     for i in range(1, len(c_array)):  # t0 时刻的数据仅用于计算t1时刻的c_hat
#         ct_array = c_array[i]
#         #print ct_array
#         ct_hat_array = c_hat_array[i - 1]
#        # print ct_hat_array
#         # ct_matrix=numpy.mat(ct_array)
#         # ct_hat_matrix=numpy.mat(ct_hat_array)
#         # delta=ct_matrix-ct_hat_matrix
#         delta = ct_array - ct_hat_array
#         d_square = numpy.square(delta)
#         d_sq_sum = numpy.sum(d_square)
#
#         data_error += d_sq_sum
#     return data_error / m
#
#
# def fitnessFuncton(c_array, lambda_sigmoid, w_array):
#     c_hat_array = getCHat(c_array, lambda_sigmoid, w_array)
#     print c_hat_array
#     data_err = objectiveFunction(c_array, c_hat_array)
#     return data_err
# end of initial version of fitness function
#####################################################
#


# selcet lines need to annotate use the short cut: ctrl+/
def initialWeightMatrix_2DArray(nGenes_int):  # generate a random weight matrix
    ini_Weight_2DArray = numpy.random.random((nGenes_int, nGenes_int))
    ini_Weight_2DArray = ini_Weight_2DArray - 0.5
    return ini_Weight_2DArray

# def main():
#     # 6 Genes(rows) , 4 time points(cols)
#     # c_list=[[0.5,0.7,0.2,0.3,0.6,0.1],[0.8,0.6,0,0.9,0.5,0.2],[0.7,0.6,0.3,0.8,0.2,0],[0.6,0.9,0.8,0,0.2,0]]
#     c_list = [[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2], [0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
#               [0.4, 0.4, 0.4, 0.4, 0.4, 0.4]]
#
#     c_array = numpy.array(c_list)
#     # 6*6 w_list
#     w_list = [[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6],
#               [1, 2, 3, 4, 5, 6]]
#     w_array = numpy.array(w_list)
#     #print fitnessFuncton(c_array, 1, w_array)
#     print fitness2(c_array, w_array,1)
#
#
# if __name__ == "__main__":
#     main()
