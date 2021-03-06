import matplotlib.pyplot as plt
import vegas
import math
import random
import graphviz
from networkx import *
from numpy import *
from celluloid import Camera


#Preparing the graph
#Num_HCW HCW nodes, each with 8 patient nodes as maximum number
Num_HCW=10
G=fast_gnp_random_graph(Num_HCW, 0.5, seed=None, directed=False)
Role=['HCW']
set_node_attributes(G, Role, 'Role')
fixed_pos=spring_layout(G)     #Setting the layout to fix the position of the HCW nodes in the plot

for n in range(len(G.nodes())):           #Adding the patient nodes, 8 for
    for i in range(8):                    #each HCW
        new=9*(n+1)+i+1
        G.add_node(new,Role=['Patient'])
        G.add_edge(n,new)

State=['Uncolonized']
set_node_attributes(G, State,'State')   #Setting all the states to uncolonized
pos = spring_layout(G,pos=fixed_pos,fixed=fixed_pos.keys())   #Layout of the node with fixed the HCW nodes


infected_list=[]  #Inizialization of the infected list
#Parameters of the model
mu_c=0.05     #Colonized discharge rate
mu_u=0.1      #Uncolonized discharge rate
alpha=0.21    #Probability of colonization
phi=0.05      #Probability that the patient added is colonized
lambd=1       #Probability to add a patient
mu_HC=0.5     #Probability of a colonized HCW to become uncolonized
p=0
#Colonized admition rate phi*lambd
#Uncolonized admition rate (1-phi)*lambd

Time=1000 #Time steps of the simulation (days)
Frac_Col=zeros(Time,'double')
Frac_Uncol=zeros(Time,'double')
Uncol_HCW=zeros(Time,'double')
Col_HCW=zeros(Time,'double')
Uncol_pat=zeros(Time,'double')
Col_pat=zeros(Time,'double')
tim=zeros(Time,'int')
#Begin of the time simulation
for t in range(Time):

    #Adding the methods to remove the diffusion after the beginning of the contagion
    if t>Time/2:
        phi=0        #Isolation of admitted colonized patient
        alpha=0.21*(1-p)    #Head washing to decrease the probability to became colonized after the contact


    #Decontaminatio of the HCW
    for i in range(Num_HCW):
        if G.nodes[i]['State']==['Colonized']:
            rand_head=random.uniform(0,1)
            if rand_head<=mu_HC:
                G.nodes[i]['State']=['Uncolonized']






    #Removing the nodes
    to_remove_list=[]
    for node in G.nodes():
        if G.nodes[node]['Role']==['Patient']:
            if G.nodes[node]['State']==['Colonized']:
                rand_colon=random.uniform(0,1)
                if rand_colon<=mu_c:
                    to_remove_list.append(node)
            if G.nodes[node]['State']==['Uncolonized']:
                rand_uncolon=random.uniform(0,1)
                if rand_uncolon<=mu_u:
                    to_remove_list.append(node)
    for j in range(len(to_remove_list)):
        G.remove_node(to_remove_list[j])

    big_node=amax(G.nodes())+1  #Setting the bigger node to begin the adding process

    #Adding the nodes kepping the maximum number of patient for each node equal to 8
    to_add_list_col=[]
    to_add_list_uncol=[]
    to_add_list_link=[]
    newnum=big_node
    for i in range(Num_HCW):                          #Cicling only on the HCW because we can only add patient to the HCW
        neighbours=list(G.neighbors(i))
        patient=[j for j in range(len(neighbours)) if G.nodes[neighbours[j]]['Role']==['Patient']] #Counting the number of patient for each node
        patient_num=len(patient)

        if patient_num<7:
            rand_admis=random.uniform(0,1)
            if rand_admis<=lambd:
                rand_pat=random.uniform(0,1)
                if rand_pat<=phi:
                    to_add_list_col.append(newnum)   #add colonized node
                    to_add_list_link.append((i,newnum))
                    newnum=newnum+1
                else:
                    to_add_list_uncol.append(newnum)  #add uncolonized node
                    to_add_list_link.append((i,newnum))
                    newnum=newnum+1



    #Adding of the nodes to the graph
    for i in range(len(to_add_list_col)):
        G.add_node(to_add_list_col[i], Role=['Patient'], State=['Colonized'])

    for i in range(len(to_add_list_uncol)):
        G.add_node(to_add_list_uncol[i], Role=['Patient'], State=['Uncolonized'])

    for i in range(len(to_add_list_link)):
        G.add_edge(*to_add_list_link[i])

    #Stop condition if the graph is empty
    if len(G.nodes())==0:
        break

    #Infection for each node
    to_infect_list=[]
    for node in G.nodes():
        if G.nodes[node]['State']==['Colonized']:
            neighbours=list(G.neighbors(node))
            for j in range(len(neighbours)):
                if G.nodes[neighbours[j]]['State']==['Uncolonized']:
                    rand_inf=random.uniform(0,1)
                    if rand_inf<=alpha:
                        to_infect_list.append(neighbours[j])


    to_infect_list=list(set(to_infect_list))  #Keeping only the unique value for each colonized node

    for i in range(len(to_infect_list)):
        G.nodes[to_infect_list[i]]['State']=['Colonized']
        infected_list.append(to_infect_list[i])

    #Counting the nodes
    for node in G.nodes():
        if G.nodes[node]['State']==['Colonized']:
            Frac_Col[t]=Frac_Col[t]+1
            if G.nodes[node]['Role']==['Patient']:
                Col_pat[t]=Col_pat[t]+1
            else:
                Col_HCW[t]=Col_HCW[t]+1
        else:
            if G.nodes[node]['Role']==['Patient']:
                Uncol_pat[t]=Uncol_pat[t]+1
            else:
                Uncol_HCW[t]=Uncol_HCW[t]+1

    Frac_Col[t]=Frac_Col[t]/len(G.nodes())
    Frac_Uncol[t]=1-Frac_Col[t]
    tim[t]=t



#Plot the ratio of colonized nodes for time step
plt.plot(tim,Frac_Col,'r',label='Colonized')
plt.plot(tim,Frac_Uncol,'g',label='Uncolonized')
plt.legend(loc='upper right')
plt.xlabel('Time')
plt.ylabel('Relative number of nodes')
plt.title('Diffusion of CPKP')
plt.show()

#Plot the total number of colonized nodes for time step
plt.plot(tim,Col_pat,'r',label='Colonized Patient')
plt.plot(tim,Uncol_pat,'g',label='Uncolonized Patient')
plt.plot(tim,Col_HCW,'y',label='Colonized HCW')
plt.plot(tim,Uncol_HCW,'blue',label='Uncolonized HCW')
plt.legend(loc='upper right')
plt.xlabel('Time')
plt.ylabel('Absolute number of nodes')
plt.title('Diffusion of CPKP')
plt.show()

