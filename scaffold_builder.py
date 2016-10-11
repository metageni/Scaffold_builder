# Scaffold_builder version 2.2

# (c)            Silva GG, Dutilh BE, Matthews TD, Elkins K, Schmieder R, Dinsdale EA, Edwards RA
# Please cite:   "Combining de novo and reference-guided assembly with Scaffold_builder", Source Code for Biology and Medicine 2013
# crAss website: http://edwards.sdsu.edu/scaffold_builder/

import string,os,sys

##########################
#  Program Defaults
##########################

parameters={"-a":95,"begin_end":0,"-p":"Scaffold","-i":80,
            "-g":5000,"-b":0,"-r":"","-q":"",
            "ambiguous":[],"subSet":[],"rawMapping":"","gaps":[],"overlapOver":0,"overlapBelow":0,"-t":300}

# [-b] behavior for rearrangements: 0 =place end-to-end;  1 = create new scaffold sequence;
# [-t] length of terminus that will be aligned
# [-i] minimum identity for merging contigs
# [-a] if this percentage of the contig length is mapped in more than one location, the contig is removed as "ambiguously mapped"
# [-g] maximum gap length allowed (default 5000nt)


#Returns the path of a program
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

##########################
# Needleman Wunsch Alignment
##########################

def needlemanWunsch(seqA,seqB,minimumIdentity):
    def ambiguityCode(sequence1,sequence2):
        iupac={'AC': 'M','AG': 'R','AT': 'W','CG': 'S','CT': 'Y','GT': 'K','AN': 'N','GN': 'N','CN': 'N','TN': 'N','NT': 'N',"TY":"C","CY":"C"}
        consensus=""
        for i in range(len(sequence1)):
            l=[sequence1[i],sequence2[i]]
            if "-" in l:l.remove("-")
            l=list(set(l))
            l.sort()

            if len(l)>1:
                letter="".join(l)
                if letter not in iupac:
                    consensus+="N"
                else:
                    consensus+=iupac[letter]
            else:
                consensus+=l[0]
        return consensus

    def createMatrix():
        matrix=[]
        for i in range(len(seqB)+1):
            matrix.append([0]*(len(seqA)+1))
        return matrix
    
    matrix=createMatrix()

    def scoreMatrix(matrix):
        column=1;maximum=[]
        for i in range(len(seqA)):
            line=1
            for j in range(len(seqB)):
                if seqB[line-1]==seqA[column-1]:
                    maximum.append(matrix[line-1][column-1]+1)
                    maximum.append(matrix[line][column-1])
                    maximum.append(matrix[line-1][column])
                else:
                    maximum.append(matrix[line-1][column-1])
                    maximum.append(matrix[line][column-1])
                    maximum.append(matrix[line-1][column])
                    
                matrix[line][column]=max(maximum)

                maximum=[]
                
                line+=1
                
            column+=1
        return matrix

    matrix=scoreMatrix(matrix)


    def perIdentity(finalA,finalB,per=0):
        for k in range(len(finalA)):
            if finalA[k] == finalB[k]:
                per+=1
                
        return  per*100.0/len(finalA) if len(finalA)>0 else 0

    def traceback(matrix,seqA,seqB):
        c=len(seqA)
        l=len(seqB)

        alignmentA="";alignmentB=""

        while (l!=0) and (c!=0):
            current=matrix[l][c]
            up=matrix[l-1][c]
            front=matrix[l][c-1]
            diagonal=matrix[l-1][c-1]

            maximo=[diagonal,front,up]

            m=max(maximo)

            i=maximo.index(m)

            if (len(set(maximo))==1) or (i==0):#diagonal
                l-=1;c-=1
                alignmentA=seqA[c]+alignmentA
                alignmentB=seqB[l]+alignmentB

            elif i==1:#front
                c-=1
                alignmentA=seqA[c]+alignmentA
                alignmentB="-"+alignmentB

            else:#i==0 #up
                l-=1
                alignmentA="-"+alignmentA
                alignmentB=seqB[l]+alignmentB
        perc=perIdentity(alignmentA,alignmentB)

        consensusSequence=""
        if perc>=minimumIdentity:consensusSequence=ambiguityCode(alignmentA,alignmentB)
            
        return [alignmentA,alignmentB,perc,consensusSequence]


    return traceback(matrix,seqA,seqB)


##########################
#  Function Definitions
##########################

#Add the fasta to a hash
def fasta2hash(fasta,my_hash={}):
    for line in "".join(file(fasta).readlines()).split(">")[1:]:
        line=line.split("\n")
        my_hash[line[0].split()[0].replace("|","")]="".join(line[1:]).upper()
    return my_hash

# Return the reverse complement of a DNA sequence
def reverseComplement(sequence):
    complement = string.maketrans('ATCGN', 'TAGCN')
    return sequence.upper().translate(complement)[::-1]

#Reverse the sequence in the hash to sequences that mapped reversed
def reverse(mapping):
    for i in mapping:
        if (("_begin" not in i[0]) and ("_end" not in i[0])) and (i[1][2]>i[1][3]):
            hashSequences[i[0]]=reverseComplement(hashSequences[i[0]])            

# To accommodate circular reference chromosomes, this function tests if the contig mapped to the start and end of the chromosome reference is identical
def beginEnd(coordHash):
    mappingBeginEnd=coordHash[parameters["begin_end"]]
    beginEnd=parameters["begin_end"]

    #BEGIN
    coordHash[beginEnd+"_begin"]=[mappingBeginEnd[0]]
    
    if (mappingBeginEnd[0][2]>mappingBeginEnd[0][3]):
        hashSequences[beginEnd+"_begin"]=reverseComplement(hashSequences[beginEnd])[mappingBeginEnd[0][3]-1:mappingBeginEnd[0][2]]
    else:
        hashSequences[beginEnd+"_begin"]=hashSequences[beginEnd][mappingBeginEnd[0][2]-1:mappingBeginEnd[0][3]]
    #END
    coordHash[beginEnd+"_end"]=[mappingBeginEnd[-1]]

    if (mappingBeginEnd[-1][2]>mappingBeginEnd[-1][3]):
        hashSequences[beginEnd+"_end"]=reverseComplement(hashSequences[beginEnd])[mappingBeginEnd[0][3]-1:mappingBeginEnd[0][2]]
    else:
        hashSequences[beginEnd+"_end"]=hashSequences[beginEnd][mappingBeginEnd[0][2]-1:mappingBeginEnd[0][3]]
    
    del coordHash[beginEnd]

    return coordHash

#Add the mummer result into a Hash Table
def coord2hash(coords):
    coords=open(coords)
    coordHash={}
    for line in coords:
        line=line.replace("|","").split()
        #Result Important Line
        if len(line)==9:
            resultInt=[int(x) for x in line[:4]]
            #is it the begin_end?
            if resultInt[0]==1:
                parameters["begin_end"]=line[-1]
            #If the sequence ID is not in the hash 
            if line[-1] not in coordHash:
                coordHash[line[-1]]=[resultInt]
            else:
                coordHash[line[-1]]+=[resultInt]       
    coords.close()

    #Seeting begin and end
    if parameters["begin_end"]!=0:coordHash=beginEnd(coordHash)

    #Store into the hash a raw copy of the mapping
    parameters["rawMapping"]=dict(coordHash)

    return coordHash

#Cleans the coords just leaving the important information to be scaffolded
# Ambiguously mapped contigs are removed:
# - contigs mapped over >X% of their length to more than one location on the reference
# - contigs mapped to more than one location on the reference on opposite strands
def cleanCoords(coordHash):
    ambiguous=[]
    #Choose the 2 best mapping to each sequence
    def chooseTop2(listResults):
        lengthList=[]
        for result in listResults:
            lengthList.append((abs(result[2]-result[3])+1,result))
        lengthList.sort()

        lengthList=lengthList[::-1][:2]
        lengthList[0]=lengthList[0][1];lengthList[1]=lengthList[1][1]

        return lengthList
        
    for sequenceID in coordHash:
        if len(coordHash[sequenceID])>1:
            #Circular // Begin-End
            if sequenceID!=parameters["begin_end"]:
                top2=chooseTop2(coordHash[sequenceID]);secondBest=top2[1]
                #AmbiguousTest  #SecondMappingLength
                if (float(abs(secondBest[2]-secondBest[3])+1)/len(hashSequences[sequenceID]))*100<parameters["-a"]:
                    coordHash[sequenceID]=top2[0]
                else:
                    parameters["ambiguous"]+=[sequenceID]
                    ambiguous.append(sequenceID)
                
        else:
            coordHash[sequenceID]=coordHash[sequenceID][0]
    #Remove from the hash the mapping which were Ambiguous
    for i in ambiguous:
        del coordHash[i]
        
    return coordHash

#Writes sequences into the final Scaffold Fasta
def writeFasta(sequenceId,sequence):
    fasta=open(parameters["-p"]+"_Scaffold.fasta","a+")
    fasta.write(">"+sequenceId+"\n"+sequence.replace("\n","")+"\n")
    fasta.close()

#Return a list with the contigs that were not mapped against the reference and write them on the final fasta file
def notMapped():
    notmapped=[]
    for i in hashSequences:
        if i not in parameters["rawMapping"] and i!=parameters["begin_end"]:
            notmapped.append(i)
            writeFasta(i+"_not_mapped",hashSequences[i])
    return notmapped

#Write the overlap and aligments
def writeOverlapAlignment(sequenceId1,sequence1,sequenceId2,sequence2,action="Overlap_"):
    overlap=open(parameters["-p"]+"_overlap_alignment/"+action+sequenceId1+"_"+sequenceId2+".fasta","w+")
    overlap.write(">"+sequenceId1+"\n"+sequence1+"\n>"+sequenceId2+"\n"+sequence2+"\n")
    overlap.close()
    
#Extend the mapping when necessary
def extendMapping(coordsList):
    for currentElement in coordsList:
        if ("_begin" not in currentElement) and ("_end" not in currentElement):
            currentLength=len(hashSequences[currentElement]);currentMapping=coordsList[currentElement]
            if (currentMapping[2]+currentMapping[3]-1)!=currentLength:
                
                if (currentMapping[2]<currentMapping[3]):
                    currentMapping[0]-=(currentMapping[2]-1)
                    if currentMapping[0]<1:currentMapping[0]=1
                    currentMapping[1]+=(currentLength-currentMapping[3])
                    
                else:
                    currentMapping[1]+=(currentMapping[3]-1)
                    currentMapping[0]-=(currentLength-currentMapping[2])
                    if currentMapping[0]<1:currentMapping[0]=1

                coordsList[currentElement]=currentMapping

    return coordsList

#Output statistics
def statistics():
    def lenght(fasta,c=0,sequences=0,longest=0,allLengths=[]):
        for line in "".join(file(fasta).readlines()).split(">")[1:]:
            line=line.split("\n")
            currentLength=len("".join(line[1:]));c+=currentLength
            allLengths.append(currentLength)
            if currentLength>longest:
                longest=currentLength
            sequences+=1
        return [c,sequences,longest,allLengths]

    def N50(Lenghts):
        #N50
        finalListLength=[]
        for i in list(set(Lenghts)):finalListLength+=[i]*Lenghts.count(i)*i
        finalListLength.sort()

        if (len(finalListLength)%2)==0:median=(finalListLength[int(len(finalListLength)/2)-1]+finalListLength[int(len(finalListLength)/2)])/2
        else:median=finalListLength[len(finalListLength)/2]
        return median
        
    #assembly length
    assembly=lenght(parameters["-q"])
    assemblyLength=assembly[0];assemblyCount=assembly[1]
    assemblyLongest=assembly[2]
    assemblyN50=N50(assembly[3])    
    #scaffold length
    scaffold=lenght(parameters["-p"]+"_Scaffold.fasta")
    scaffoldLength=scaffold[0];scaffoldCount=scaffold[1]
    scaffoldLongest=scaffold[2]
    scaffoldN50=N50(scaffold[3])

    return [assemblyLength,scaffoldLength,assemblyCount,scaffoldCount,assemblyLongest,scaffoldLongest,assemblyN50,scaffoldN50]
    
#Write the output with the scaffold information
def writeOutput(finalScaffold):
    info="Scaffold_builder version 2.2 log file\n\nReference sequence: "+parameters["-r"]+"\n\n"

    #beginEnd?
    if parameters["begin_end"]!=0:info+=parameters["begin_end"]+" maps to the beginning and end of the reference, suggesting that this reference is circular.\n\n"

    #Contigs ambiguous
    if len(parameters["ambiguous"])>1:
        info+="Contigs mapped ambiguously:"+" The contigs below map to the reference at least 2 times over  >= "+str(parameters["-a"])+"% of it's length, so it is not scaffolded.\n"
        for i in parameters["ambiguous"]:
            writeFasta(i+"_ambiguous_mapping",hashSequences[i])
            info+=i+"\n"  
            
    #SubSet
    if len(parameters["subSet"])>1:
        info+="\nContigs mapped entirely within a region where another contig has already been scaffolded are removed:\n\n"
        for i in parameters["subSet"]:
            temp=[]
            for n in i:
                temp.append(n[0]+" ("+str(n[1][0])+"-"+str(n[1][1])+")")

            subset=temp[0].split(" (")[0];writeFasta(subset+"_overlapping_hits_sub_region",hashSequences[subset])
            info+=temp[0]+" lies within the region of "+temp[1]+"\n"

    #Writing Final Scaffold
    info+="\nFinal scaffolding:\n\n"
    info+="#contig\t\t\t\t\t5'ref\t\t\t\t\t3'ref\t\t\t\t\t5'contig\t\t\t\t\t3'contig\t\t\t\t\tlength\n"
    for i in finalScaffold:
        info+=str("\t\t\t\t\t".join([i[0]]+[str(x) for x in i[1]]+[str(abs(i[1][-1]-i[1][-2])+1)]))+"\n"

    #NotMapped
    info+="\nContigs that were not mapped to the reference by Nucmer\n"
    for i in notMapped():
        info+=i+"\n"
    
    statistic=statistics()
    
    #Statistics
    info+="\nScaffolding statistics:\n\t\t\t\t\tAssembly\tScaffold\n"
    info+="Total length\t"+str(statistic[0])+"\t"+str(statistic[1])+"\n"
    info+="Number of sequences\t"+str(statistic[2])+"\t"+str(statistic[3])+"\n"
    info+="Average length\t"+str(round(statistic[0]/statistic[2],2))+"\t"+str(round(statistic[1]/statistic[3],2))+"\n"
    info+="Longest contig\t"+str(statistic[4])+"\t"+str(statistic[5])+"\n"
    info+="N50\t"+str(statistic[6])+"\t"+str(statistic[7])+"\n"
    info+="Non-overlapping contig pairs\t"+str("-")+"\t"+str(len(parameters["gaps"]))+"\n"
    info+="Total length of gaps\t"+str("-")+"\t"+str(sum(parameters["gaps"]))+"\n"
    if len(parameters["gaps"])>0:
        info+="Longest gap\t"+str("-")+"\t"+str(max(parameters["gaps"]))+"\n"
        info+="Shortest gap\t"+str("-")+"\t"+str(min(parameters["gaps"]))+"\n"
        info+="Average gap size\t"+str("-")+"\t"+str(round(float(sum(parameters["gaps"]))/len(parameters["gaps"]),2))+"\n"
    else:
        info+="Longest gap\t"+str("-")+"\t0\n"
        info+="Shortest gap\t"+str("-")+"\t0\n"
        info+="Average gap size\t"+str("-")+"\t0\n"

    info+="Overlapping contig pairs (>=80% id)\t-\t"+str(parameters["overlapOver"])+"\n"
    info+="Overlapping contig pairs (<80% id)\t-\t"+str(parameters["overlapBelow"])

    #Write in the output everything
    output=open(parameters["-p"]+"_output.txt","w+")
    output.write(info)
    output.close()
    
#Scaffold the clean mapping
def scaffold(mapping):
    c=0;scaffold="";scaffoldNumber=1
    while c<(len(mapping)-1):
        currentElement=mapping[c]
        nextElement=mapping[c+1]

        #Gap Length | If it is negative it is a overlap, otherwise, GAP
        checkGap=(nextElement[1][0]-currentElement[1][1]-1)
        
        if checkGap>=0:#Check how big is the gap and if we should insert Ns
            #Store the length to use as statistics later
            parameters["gaps"]+=[checkGap]
            #Break 
            if checkGap>=parameters["-g"] and len(scaffold)>0:
                writeFasta("Scaffold_"+str(scaffoldNumber),scaffold);scaffold=""
                scaffoldNumber+=1

            else:
                # Filling gaps between non-overlapping adjacently mapped contigs with the appropriate number of Ns
                scaffold+=hashSequences[currentElement[0]]+"N"*checkGap
        else:# If the contig is not the last one mapped to the reference and the next one overlaps with it, they should be aligned
            overlapLength=abs(checkGap)

          
            if overlapLength>parameters["-t"]:overlapLength=parameters["-t"]
            # Extract the termini of the adjacently mapped contigs
            currentOverlap=hashSequences[currentElement[0]][-overlapLength:]
            nextOverlap=hashSequences[nextElement[0]][:overlapLength]
            
            currentNoOverlap=hashSequences[currentElement[0]][:-overlapLength]
            nextNoOverlap=hashSequences[nextElement[0]][overlapLength:]

            # Align using the needlemanWunsch Algorithm
            alignment=needlemanWunsch(currentOverlap,nextOverlap,parameters["-i"])
            
            #If the aligment didnt have >= to (alignment[-2]<parameters["-i"]) we try to find it 
            if (alignment[-2]<parameters["-i"]) and (overlapLength<parameters["-t"]) and (overlapLength+10)<(parameters["-t"]+1):
                results=[]
                for i in range(overlapLength+10,parameters["-t"]+1,10):
                    currentOverlap=hashSequences[currentElement[0]][-i:]
                    nextOverlap=hashSequences[nextElement[0]][:i]
                    
                    newAlignment=needlemanWunsch(currentOverlap,nextOverlap,parameters["-i"])

                    results.append((newAlignment[-2],newAlignment,i))

                results.sort()

                results=results[::-1][0]
                if results[0]>alignment[-2]:
                    
                    alignment=results[1];overlapLength=results[2]
                    # Update the variables if we found a better length with a better Identitify
                    currentOverlap=hashSequences[currentElement[0]][-overlapLength:]
                    nextOverlap=hashSequences[nextElement[0]][:overlapLength]
                    
                    currentNoOverlap=hashSequences[currentElement[0]][:-overlapLength]
                    nextNoOverlap=hashSequences[nextElement[0]][overlapLength:]

            #Get the consensus sequence || currentSequence(w/o overlap part) + consensus + Update hash with nextSequence(w/o overlap part)
            if alignment[2]>=parameters["-i"]:
                #statistics
                parameters["overlapOver"]+=1
                scaffold+=currentNoOverlap+alignment[3]
                hashSequences[nextElement[0]]=nextNoOverlap
            #Add the current Sequence into the scaffold
            else:
                #statistics
                parameters["overlapBelow"]+=1
                scaffold+=hashSequences[currentElement[0]]
                
                ##Break the scaffold if the identity is less than 80% (default)
                if parameters["-b"]!=0 and len(scaffold)>0:
                    writeFasta("Scaffold_"+str(scaffoldNumber),scaffold);scaffold=""
                    scaffoldNumber+=1
                
            #REPORT the overlaps to the user
            writeOverlapAlignment(currentElement[0],currentOverlap,nextElement[0],nextOverlap)
            writeOverlapAlignment(currentElement[0],alignment[0],nextElement[0],alignment[1],"Alignment_")
            
        c+=1
    try:
        scaffold+=hashSequences[nextElement[0]]
    except:
        pass
    
    if len(scaffold)>0:
        writeFasta("Scaffold_"+str(scaffoldNumber),scaffold)

    #Writing Output
    writeOutput(mapping)
    
#Sort the scaffold aka set it to be scaffolded
def sortMapping(coordHash):
    #Remove the subSet mapping
    def removeSubSet(sortedList):
        c=0;subSet=[]
        while c<(len(sortedList)-1):
            if ((sortedList[c][1][0])<(sortedList[c+1][1][0])) and ((sortedList[c][1][1])>(sortedList[c+1][1][1])):
                parameters["subSet"]+=[[sortedList[c+1],sortedList[c]]]
                sortedList.remove(sortedList[c+1])
                c-=1
            c+=1
        return sortedList


    #Extend Mapping
    coordHash=extendMapping(coordHash)
        
    
    mappingSorted=[]
    for element in coordHash:
        mappingSorted.append((coordHash[element][0],[element,coordHash[element]]))
    mappingSorted.sort()
    
    t=mappingSorted;mappingSorted=[]

    for i in t:
        mappingSorted.append(i[1])

    #Reverse sequences thast neeed it
    reverse(mappingSorted)
    
    #Remove SubSet
    mappingSorted=removeSubSet(mappingSorted)

    #### SCAFFOLD
    scaffold(mappingSorted)


helpMsg="\nScaffold_builder version v 2.2\n\nUsage:\npython scaffold_builder.py -q query_contigs.fna -r reference_genome.fna -p output_prefix [-t] [-i] [-a] [-b]\n\n\t-q fasta file of contigs\n\t\tRequired. Query contigs in Fasta format. These contigs may be the output of a de novo\n\t\tassembly program such as Newbler, Velvet or MIRA.\n\n\t-r fasta file containing reference genome\n\t\tRequired. Reference genome in Fasta format. This should preferably be a completed genome\n\t\tsequence.\n\n\t-p prefix output files\n\t\tRequired. All the output files have this project name as prefix.\n\n\t-t length of terminus that will be aligned (default "+str(parameters["-t"])+" nt)\n\t\tAt any break between two contigs, scaffold_builder checks whether the termini\n\t\tof the adjacent contigs are homologous by aligning them using Smith-Waterman's Algorithm, and\n\t\tcombines them if that is the case.\n\n\t-i minimum identity for merging contigs (default "+str(parameters["-i"])+"%)\n\t\tIf the termini are similar, scaffold_builder assumes that the contigs should\n\t\thave been combined by the assembly program, but the similarity was probably\n\t\tbelow the assembly thresholds, or the contigs were not merged due to ambiguous\n\t\tread mapping. The sequences are combined and in the case that non-identical\n\t\tnucleotides are aligned, the IUPAC code of their consensus is placed in the\n\t\tresulting sequence.\n\n\t-a minimum length for ambiguously mapped contigs (default "+str(parameters["-a"])+"%)\n\t\tIf a contig maps to more than one location on the reference genome, it will\n\t\tnot be scaffolded because it's location is ambiguous. This parameter defines\n\t\thow much of the length of a contig should be mapped in more than one location\n\t\tfor it to be considered ambiguously mapped.\n\n\t-b 0/1 dictates behavior for rearrangements (default "+str(parameters["-b"])+")\n\t\t0: place end-to-end\n\t\t1: create new scaffold sequence\n\t\tIf the mapping of the contigs onto the reference suggests that they overlap,\n\t\tbut the contig termini are too dissimilar to join them, this option dictates\n\t\twhether scaffold_builder places the contigs end-to-end (default; deletions\n\t\texpected) or to start a new scaffold sequence (inversions expected).\n\n\t-g maximum gap length allowed (default "+str(parameters["-g"])+"nt)"

def run():
    #Set parameters
    userParameters=sys.argv[1:]
    if "-h" in userParameters:
        print userParameters
        return 0
    # Check if the user missed any essential parameter
    if ("-q" not in userParameters):
        print "#Please inform the query file -q"

    if ("-r" not in userParameters):
        print "#Please inform the reference file -r"
    else:
        for i in range(0,len(userParameters),2):
            parameters[userParameters[i]]=userParameters[i+1]
        parameters["-t"]=int(parameters["-t"])
        parameters["-g"]=int(parameters["-g"])
        parameters["-i"]=int(parameters["-i"])
        parameters["-a"]=int(parameters["-a"])
        parameters["-b"]=int(parameters["-b"])
        
        #Run Nucmer // Creates the folder to the aligments and overlaps
        os.system(which("nucmer")+" "+parameters["-r"]+" "+parameters["-q"]+" && "+which("show-coords")+" out.delta >"+parameters["-p"]+".coords && rm out.delta && mkdir "+parameters["-p"]+"_overlap_alignment")
        return True


try:
    if run()==True:
        #Main
        hashSequences=fasta2hash(parameters["-q"])
        coords=coord2hash(parameters["-p"]+".coords")
        cleanCoords=cleanCoords(coords)
        sortMapping(cleanCoords)

        print "Scaffolded :)"
except:
    print helpMsg
