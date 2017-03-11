import random
import numpy as np
import math
import os
from scipy import  ndimage
def continus_calculate(b):
    chicun = int(math.sqrt(len(b)))
    b[pow(chicun, 2) - 1] = 1
    b[pow(chicun, 2) - 1] = 1
    b[pow(chicun, 2) - chicun - 1] = 1
    p = [x + 1 for x in range(len(b)) if b[x] == 1]
    position_matrix = np.zeros([chicun + 2, chicun + 2])
    for x in range(len(p)):
        position_matrix[int((p[x] - 0.001) / chicun) + 1, int((p[x] - 0.001) % chicun) + 1] = 1
    for i in range(chicun-1,0,-1):
        for j in range(chicun,i-1,-1):
            if not(i==j):
                # print i,j
                if position_matrix[i,j]==1:
                    if  (position_matrix[i,j+1]==0 and position_matrix[i+1,j+1]==0 and position_matrix[i+1,j]==0):
                        position_matrix[i, j]=0
                if position_matrix[j,i]==1:
                    if  (position_matrix[j,i+1]==0 and position_matrix[j+1,i+1]==0 and position_matrix[j+1,i]==0 and position_matrix[j-1,i+1]==0):
                        position_matrix[j, i]=0
            if i==j:
                if position_matrix[j,i]==1:
                    if  (position_matrix[j,i+1]==0 and position_matrix[j+1,i+1]==0 and position_matrix[j+1,i]==0):
                        position_matrix[j, i]=0
    # print(position_matrix)
    return (sum(sum(position_matrix)))
def continues_format(b,volume,iter):
    chicun=int(math.sqrt(len(b)))
    b[pow(chicun,2)-1]=1
    b[pow(chicun,2)-1]=1
    b[pow(chicun,2)-chicun-1]=1
    p=[x+1 for x in range(len(b)) if b[x]==1]
    position_matrix=np.zeros([chicun+2,chicun+2])
    for x in range(len(p)):
        position_matrix[int((p[x]-0.001)/chicun)+1,int((p[x]-0.001)%chicun)+1]=1
    position_matrix_copy=position_matrix.copy()
    # position_matrix=np.array([[0,0,0,0,0,0,0],[0,1,1,1,0,0,0],[0,1,0,1,0,0,0],[0,1,0,1,0,0,0], \
    #     [0, 1, 1, 1, 1, 0, 0],[0,0,0,0,1,1,0],[0,0,0,0,0,0,0]])


    # print('===========position_matrix_copy')
    # print(position_matrix_copy)
    # print(sum(sum(position_matrix_copy)))

    for i in range(chicun-1,0,-1):
        for j in range(chicun,i-1,-1):
            if not(i==j):
                # print i,j
                if position_matrix[i,j]==1:
                    if  (position_matrix[i,j+1]==0 and position_matrix[i+1,j+1]==0 and position_matrix[i+1,j]==0):
                        position_matrix[i, j]=0
                if position_matrix[j,i]==1:
                    if  (position_matrix[j,i+1]==0 and position_matrix[j+1,i+1]==0 and position_matrix[j+1,i]==0 and position_matrix[j-1,i+1]==0):
                        position_matrix[j, i]=0
            if i==j:
                if position_matrix[j,i]==1:
                    if  (position_matrix[j,i+1]==0 and position_matrix[j+1,i+1]==0 and position_matrix[j+1,i]==0):
                        position_matrix[j, i]=0

    # print('============position_matrix')
    # print(position_matrix)
    # print(sum(sum(position_matrix)))


    # print('============position_matrix_copy-position_matrix')
    # print(position_matrix_copy-position_matrix)
    # print(sum(sum(position_matrix_copy-position_matrix)))


    position_matrix[len(position_matrix) - 1, :] = 1
    position_matrix[:, len(position_matrix) - 1] = 1

    position_matrix_fill=ndimage.binary_fill_holes(position_matrix).astype(int)

    position_matrix_fill[len(position_matrix_fill) - 1, :] = 0
    position_matrix_fill[:, len(position_matrix_fill) - 1] = 0
    fill_sum = sum(sum(position_matrix_fill))

    # print('================position_matrix_fill===befor volume constrain')
    # print(position_matrix_fill)
    # print(sum(sum(position_matrix_fill)))

    for i in range(iter):
        if fill_sum>volume:
            position_matrix_fill[len(position_matrix_fill)-1,:]=1
            position_matrix_fill[:,len(position_matrix_fill)-1] = 1
            contin_vertor=[]
            for i in range(1,len(position_matrix_fill)-1):
                for j in range(1,len(position_matrix_fill)-1):
                    if position_matrix_fill[i,j]==1:
                        t=[]
                        count_border_number = 0
                        for ii in range(-1,2,1):
                            for jj in range(-1,2,1):
                                if position_matrix_fill[i+ii,j+jj]==0:
                                    count_border_number+=1
                        t.append(count_border_number)
                        t.append(i)
                        t.append(j)
                        contin_vertor.append(t)
            mat1 = np.array(contin_vertor)
            mat1 = mat1[mat1[:, 0].argsort()]
            for i in range(fill_sum-1,0,-1):
                continus_calculate_list_befor=[]
                for i_1 in range(1, len(position_matrix_fill) - 1):
                    for j_1 in range(1, len(position_matrix_fill) - 1):
                        continus_calculate_list_befor.append(int(position_matrix_fill[i_1, j_1]))
                c1=continus_calculate(continus_calculate_list_befor)

                if c1 < volume:
                    break
                
                position_matrix_fill[mat1[i,1],mat1[i,2]]=0

                continus_calculate_list_after=[]
                for i_2 in range(1, len(position_matrix_fill) - 1):
                    for j_2 in range(1, len(position_matrix_fill) - 1):
                        continus_calculate_list_after.append(int(position_matrix_fill[i_2, j_2]))
                c2=continus_calculate(continus_calculate_list_after)
                if c1-c2>1:
                    position_matrix_fill[mat1[i, 1], mat1[i, 2]] = 1


            position_matrix_fill[len(position_matrix_fill)-1,:]=0
            position_matrix_fill[:,len(position_matrix_fill)-1]=0
            fill_sum=sum(sum(position_matrix_fill))

    # print('================position_matrix_fill===after volume constrain')
    # print(position_matrix_fill)
    # print('--------')
    # print(sum(sum(position_matrix_fill[1:len(position_matrix_fill)-1,1:len(position_matrix_fill)-1])))

    b[:]=[]
    for i in range(1,len(position_matrix_fill)-1):
        for j in range(1,len(position_matrix_fill)-1):
            b.append(int(position_matrix_fill[i,j]))
    # print(b)


if __name__=="__main__":
    chicun=20
    b=[]
    for i in range(pow(chicun,2)):
        a=random.randint(0,1)
        b.append(a)
    print(b)
    continues_format(b,40,30)