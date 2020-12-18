from biopandas.pdb import PandasPdb
import math
pdbInput = input('Please enter the name of the .pdb file with the extension (and the corresponding directory), e.g., 2bti.pdb, ./assignment3_proteins/3r0a.pdb etc.\n')

ppdb = PandasPdb()
ppdb.read_pdb(pdbInput)
#ppdb.read_pdb('4ix3.pdb')
ddf = ppdb.df['ATOM']
dddf = ddf.drop(["record_name", "charge", "line_idx"], axis=1)
# Boolean variables to select the correct coordinates of the atoms.
alpha = (dddf['residue_name'] == "GLY") & (dddf['atom_name'] == "CA") 
beta = (dddf['atom_name'] == "CB") & (dddf['residue_name'] != "GLY")
coorDdf = dddf[alpha | beta] 

chainA = coorDdf['chain_id'] == "A" # Boolean variable for chain ids with A.
chainB = coorDdf['chain_id'] == "B" # Boolean variable for chain ids with B.
ddfA = coorDdf[chainA] # Dataframe for A
ddfB = coorDdf[chainB] # Dataframe for B

# Converting the dataframes into lists.
lddfA = ddfA.values.tolist()
lddfB = ddfB.values.tolist()
# List of tuples that will contain interacting residue pairs.
pairs = []

for atomA in lddfA:
    for atomB in lddfB:
        distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(atomA[10:13], atomB[10:13])]))
        if distance < 8:
            pairs.append((atomA,atomB))
          
pairSize = len(pairs)
pairsCop = pairs.copy()
temp = pairsCop[0]
pairsCop.remove(temp)
eachGroup = [] # List for each single group ing groups list.
eachGroup.append(temp)
groups = [] # Group list for interaction hot spots.
interactList = [] # Boolean list which consists of interacting values to compare distances.
while len(pairsCop) != 0:
   
    for i in range(0,len(pairsCop)):
        if temp != pairsCop[i]: # if temp is the last element, we can't make further computation.
            interactList = []
            for j in range(0,len(eachGroup)):            
                first = eachGroup[j][0]
                second = eachGroup[j][1]
                nextF = pairsCop[i][0]
                nextS = pairsCop[i][1]
                dist1 = math.sqrt(sum([(a - b) ** 2 for a, b in zip(first[10:13], nextF[10:13])]))
                dist2 = math.sqrt(sum([(a - b) ** 2 for a, b in zip(first[10:13], nextS[10:13])]))
                dist3 = math.sqrt(sum([(a - b) ** 2 for a, b in zip(second[10:13], nextF[10:13])]))
                dist4 = math.sqrt(sum([(a - b) ** 2 for a, b in zip(second[10:13], nextS[10:13])]))
                interacting = (dist1 < 8) or (dist2 < 8) or (dist3 < 8) or (dist4 < 8)
                interactList.append(interacting)
            if any(interactList):
                temp = pairsCop[i]
                eachGroup.append(temp)
                pairsCop.remove(temp)
                break
    
        else: # then temp is the last element in the pairsCop and it is not in the same group with the rest.
            pairsCop.remove(temp)        
    # if it is not interacting with the others -> new interaction hotspot(group) is needed.       
    if all(distBool == False for distBool in interactList) and i == (len(pairsCop)-1): 
        temp = pairsCop[0]
        groups.append(eachGroup)
        eachGroup = []
        eachGroup.append(temp)
        if len(pairsCop) != 1: # if there is MORE THAN one element in the pairsCop
            pairsCop.remove(temp)
    elif len(pairsCop) == 0:
        groups.append(eachGroup)

groupSize = len(groups)
print ("There are " + str(pairSize) + " interacting pairs.")
for i in range(0,groupSize):
    for j in range(0,len(groups[i])):
        print ("Group " + str(i+1) + ": " + str(groups[i][j][0][4]) + "(" + str(groups[i][j][0][7]) + ")-" + str(groups[i][j][1][4]) + "(" + str(groups[i][j][1][7]) + ")" )
print ("Number of groups = " + str(groupSize))    
