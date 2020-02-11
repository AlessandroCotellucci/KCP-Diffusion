import matplotlib.pyplot as plt
import vegas
import math
import random
import graphviz
from networkx import *
from numpy import *
from celluloid import Camera

#Preparation the animation
fig=plt.figure()
camera=Camera(fig)

#Preparing the graph
#Num_HCW HCW nodes, each with 8 patient nodes as maximum number
Num_HCW=10
G=fast_gnp_random_graph(Num_HCW, 0, seed=None, directed=False)
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
mu_c=0.3     #Colonized discharge rate
mu_u=0.3     #Uncolonized discharge rate
alpha=0.4    #Probability of colonization
phi=0.3      #Probability that the patient added is colonized
lambd=1      #Probability to add a patient
mu_HC=0.5    #Probability of a colonized HCW to become uncolonized
#Colonized admition rate phi*lambd
#Uncolonized admition rate (1-phi)*lambd

Time=100 #Time steps of the simulation (days)

#Begin of the time simulation
for t in range(Time):

    #plot at each step
    node_colorP=[]
    node_HCW=[]
    node_colorH=[]
    node_Patient=[]
    for node in G.nodes():                          #Cicle to set the color for each node
        if G.nodes[node]['State']==['Uncolonized']: #green for uncolonized patient and
            if G.nodes[node]['Role']==['HCW']:      #red for colonized patient.
                node_colorH.append('green')         #Also the style is setted: points for
                node_HCW.append(node)               #patient and triangles for the HCW
            if G.nodes[node]['Role']==['Patient']:
                node_colorP.append('green')
                node_Patient.append(node)
        if G.nodes[node]['State']==['Colonized']:
            if G.nodes[node]['Role']==['HCW']:
                node_HCW.append(node)
                node_colorH.append('red')
            if G.nodes[node]['Role']==['Patient']:
                node_Patient.append(node)
                node_colorP.append('red')

    pos = spring_layout(G,pos=pos,fixed=fixed_pos.keys())
    f=[]
    red,=plt.plot(f,"rs", markersize=10)
    green,=plt.plot(f,"gs",markersize=10)
    patient,=plt.plot(f,"bo",markersize=10)
    HCW,=plt.plot(f,"b^",markersize=10)
    plt.title('Diffusion KCP')
    plt.legend([red,green,patient,HCW],["Colonized","Uncolonized","Patient","HCW"])
    draw(G, with_labels=False, pos=pos,nodelist=node_HCW, node_color=node_colorH, node_shape='^')
    draw(G, with_labels=False, pos=pos,nodelist=node_Patient, node_color=node_colorP, node_shape='o')
    camera.snap()


    #Moving the HCW to the next room for the next cycle
    Node_state=get_node_attributes(G,'State')
    for i in range(Num_HCW):
        G.nodes[i]['State']=Node_state[(i+3)%Num_HCW]


    #Adding the methods to remove the diffusion after the beginning of the contagion
    if t>15:
        phi=0        #Isolation of admitted colonized patient
        alpha=0.3    #Antibiotic restriction decreasing alpha
        #Hand washing for HCW
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



#Final Plot
node_colorP=[]
node_HCW=[]
node_colorH=[]
node_Patient=[]
for node in G.nodes():
    if G.nodes[node]['State']==['Uncolonized']:
        if G.nodes[node]['Role']==['HCW']:
            node_colorH.append('green')
            node_HCW.append(node)
        if G.nodes[node]['Role']==['Patient']:
            node_colorP.append('green')
            node_Patient.append(node)
    if G.nodes[node]['State']==['Colonized']:
        if G.nodes[node]['Role']==['HCW']:
            node_HCW.append(node)
            node_colorH.append('red')
        if G.nodes[node]['Role']==['Patient']:
            node_Patient.append(node)
            node_colorP.append('red')

pos = spring_layout(G,pos=pos,fixed=fixed_pos.keys())
red,=plt.plot(f,"rs", markersize=10)
green,=plt.plot(f,"gs",markersize=10)
patient,=plt.plot(f,"bo",markersize=10)
HCW,=plt.plot(f,"b^",markersize=10)
plt.title('Diffusion KCP')
plt.legend([red,green,patient,HCW],["Colonized","Uncolonized","Patient","HCW"])
draw(G,pos,nodelist=node_HCW, node_color=node_colorH, node_shape='^')
draw(G,pos,nodelist=node_Patient, node_color=node_colorP, node_shape='o')
camera.snap()
animation=camera.animate()
animation.save(filename='KCP_enhanced_move.html')
plt.show()
