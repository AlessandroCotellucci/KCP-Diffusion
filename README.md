# CPKP Diffusion on a network
Implementation on a network of the model for the diffusion of CPKP proposed in the article (Transmission Dynamics of Carbapenemase-Producing
Klebsiella Pneumoniae and Anticipated Impact of Infection Control Strategies in a Surgical Unit) with some change.
The model implemented allow also the contagion between HCW. The model is based on a HCW cluster each HCW node has a fixed number of patient.\
Bibliography:
1. Transmission Dynamics of Carbapenemase-Producing Klebsiella Pneumoniae and Anticipated Impact of Infection Control Strategies in a Surgical Unit; Vana Sypsa, Mina Psichogiou, Georgia-Aikaterina Bouzala, Linos Hadjihannas, Angelos Hatzakis, Georgios L. Daikos
2. Mathematics of Epidemics on Networks; István Z. Kiss, Joel C. Miller, Péter L. Simon
3. Model versions and fast algorithms for network epidemiology; Petter Holme


Image description:
1. Barabasi_Albert_Diffusion.png: Relative number of contaminated nodes respect to the p (percentage of head washing in the HCW group) for the HCW's network made by a Barabasi Albert graph, with high probability admission
2. Initial_Barabasi_albert_Graph.png: Plot of the graph with Barabasi Albert structure
3. Initial_Moving_HCW_Graph.png: Plot of the graph with moving HCW graph structure (at each time step the HCWs are all shifted to the next patient room)
4. Initial_Random_Graph.png: Plot of the graph with random structure
5. Maximum_contamination_restriction_(RGP=0.7).png: Plot of the ratio of contaminated nodes and uncontaminated nodes respect to che time in the high contamination restriction (p=0.7)
6. Maximum_contamination_restriction_(RGP=0.7)_Absolute_num.png: Plot of absolute value of nodes respect to che time in the high contamination restriction (p=0.7)
7. Minimum_contaminatio_restriction_(RG).png: Plot of the ratio of contaminated nodes and uncontaminated nodes respect to che time in the low contamination restriction (p=0)
8. Minimum_contaminatio_restriction_(RG)_absolute.png: Plot of absolute value of nodes respect to che time in the low contamination restriction (p=0.7)
9. Moving_HCW_Graph_Diffusion.png: Relative number of contaminated nodes respect to the p (percentage of head washing in the HCW group) for the HCW's network made by a moving HCW graph (at each time step the HCWs are all shifted to the next patient room) graph, with high probability admission
10. Random_Graph_Diffusion.png: Relative number of contaminated nodes respect to the p (percentage of head washing in the HCW group) for the HCW's network made by a random graph, with high probability admission
11. Mean_field_T.png: mean field threshold
12. Mean_field_p=07.png: mean field contagion restriction p=0.7
