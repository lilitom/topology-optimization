import random
def Formatpoprandom(pop,ones):
    bit_1_list=[]
    for k in range(len(pop)):
        bit_1_list=[i for i in range(len(pop[k])) if pop[k][i]==1]
        if len(bit_1_list)>ones:
            sclice=random.sample(bit_1_list,len(bit_1_list)-ones)
            for j in sclice:
                pop[k][j]=0
pop=[[1,1,1,1,1,1,1,1],[0,1,0,1,1,0,1,1]]
Formatpoprandom(pop,3)
print(pop)