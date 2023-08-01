
import argparse as ap
import numpy as np

class SmithWaterman():
    # initializer method
    def __init__(self, s1, s2, m, mm, g):
        self.__seq1=s1.upper() # first sequence to align
        self.__seq2=s2.upper() # second sequence to align 
        self.__match=m # match score
        self.__mismatch=mm # miscmatch score
        self.__gap=g # gap score
        ## Initialization of two matrix -> need to be len(sequences) +1 because the forst row and the first column are all zeros -> because equal to all gaps
        self.matrix_scores=np.zeros( (len(s1)+1, len(s2)+1) ) # matrix to create the score for the process of filling with sthe best scores # nb sequence 1 in the row, seq 2 in the columns
        self.matrix_backtracking=np.chararray( (len(s1)+1, len(s2)+1), itemsize=9) # matrix to collect the possible paths from wich we obtained the score in the box cell -> needed for the process of backtracking
        for i in range(len(s1)+1): # seq 1 is the row, seq 2 is the columns
            for j in range(len(s2)+1):
                self.matrix_backtracking[i,j]=""
        self.dict_sequences={} # i create the object that i will use as container for teh alignment i will generate later on in the algorithm
        self.list_alig=[""]
        self.temp=0

    def getSeq1(self): ## get the sequence 1 that is masked!
        return(self.__seq1)
    
    def getSeq2(self): ## get the sequence 2 that is masked!
        return(self.__seq2)
    
    def getMatch(self): ## get the match score
        return(self.__match)
    
    def getMisMatch(self): ## get the mismatch score
        return(self.__mismatch)
    
    def getGap(self): ## get the gap value
        return(self.__gap)
    
    ## Matrix Filling 
    ## The second of the algorithm is filling the entire matrix, so it is more important to know the neighbor values (diagonal, upper and left) of the current cell to fill each and every cell.
    def M_Filling(self):
        sequence1=self.getSeq1()
        sequence2= self.getSeq2()
        s_match=self.getMatch()
        s_mismatch=self.getMisMatch()
        s_gap=self.getGap()
        # print(self.matrix_backtracking)
        # print(self.matrix_scores)
        for i in range(1, len(sequence1)+1): # itero sulle righe della matrix -> from row 1 till row len(seq1)+1 perchè non considero la riga 0 visto che sarebbero tutti 0 visti i gap.
            for j in range(1, len(sequence2)+1): # itero sulle colonne -> from column 1 till row len(seq2)+1
                #self.matrix_backtracking[i,j]="" # i create a string for every position of the backtracking matrix, will be needed after to insert the paths (="arrows")
                # print("coordinates \t",i, j)
                # matrix_backtarcking[i,j] returns the element at row i and column j
                upper=self.matrix_scores[i-1,j]+s_gap # score if we have a gap coming from up
                left=self.matrix_scores[i,j-1] + s_gap # score if we come from the left
                if sequence1[i-1]==sequence2[j-1]: # check if the sequences in the position i and j have thE same nucleotide therefore there is a match
                    diagonal=self.matrix_scores[i-1, j-1]+s_match
                else: # otherwise there is a mismatch
                    diagonal=self.matrix_scores[i-1, j-1]+s_mismatch
                
                self.matrix_scores[i,j]=max(diagonal, upper, left,0)
                

                ## also we create the matrix that i will be using for the process of backtracking with the paths
                list_scores= [diagonal, upper, left, 0]
                
                number_max=list_scores.count(self.matrix_scores[i,j]) # i count how many time I have the max among the possible values for teh score -> this is needed because i can obtain the same score from more than one position among upper, left or diagonal and i need to remember this!
                position_max=list_scores.index(int(self.matrix_scores[i,j]))
                while number_max>0: # when i get to 0 i need to stop because i don't have any repetition nomore
                    #print(number_max)
                     # i cancel one repetition in the number of max score present
                    # print(position_max, number_max)
                    #print(type(self.matrix_backtracking[i,j]))
                    if position_max==0: # if the max correspond to 0 this means we have a match or a mismatch, we come from the diagonal -> because diagonal is the first element of list_scores
                        if type(self.matrix_backtracking[i,j]) == str:
                            self.matrix_backtracking[i,j]= str(self.matrix_backtracking[i,j])+"D" # because we come from the diagonal
                        else:
                            self.matrix_backtracking[i,j]= self.matrix_backtracking[i,j].decode("utf-8")+"\tD" # i use decode to convert the byte class to string and then add the value for that path
                    elif position_max==1:
                        if type(self.matrix_backtracking[i,j]) == str:
                            self.matrix_backtracking[i,j]= str(self.matrix_backtracking[i,j])+"U" # we come from upper
                        else:
                            self.matrix_backtracking[i,j]= self.matrix_backtracking[i,j].decode("utf-8")+"\tU"
                    elif position_max==2:
                        if type(self.matrix_backtracking[i,j]) == str:
                            self.matrix_backtracking[i,j]= str(self.matrix_backtracking[i,j])+"L" # we come from left
                        else:
                            self.matrix_backtracking[i,j]= self.matrix_backtracking[i,j].decode("utf-8")+"\tL"
                    
                    number_max=number_max-1
                    if number_max !=0 and position_max != len(list_scores):
                        position_max=list_scores.index(int(self.matrix_scores[i,j]), position_max+1) # i find the position of the maximum in the list of the scores, this way i can understand where i am from

                    # if the max correspond to 0, since it has position 3 in list_scores, we do nothing, no backtracking, no path needed

        #print(self.matrix_scores)
        #print(self.matrix_backtracking)

    def greedy(self):
        max=np.max(self.matrix_scores) # return the max value of the score in the matrix 
        index=np.where(self.matrix_scores==max) # return the coordinates of the max in the matrix -> possible using the max value found early -> NB return an arrat object not the coordinates
        new_index=tuple(zip(*index)) # convert the array type to tuple to get faster the coordinates of the matrix in which there is a max value


        for tupletta in new_index: # iterate over the coordinates we have for teh max value
            # print(self.matrix_scores[tupletta])
            # print(self.matrix_backtracking[tupletta])
            self.temp=0
            self.list_alig=[""]           

            sw.backtracking(tupletta) # i call the backtracking here
        
        

    def backtracking(self, coord):
        #print(coord, self.matrix_scores[coord], self.matrix_backtracking[coord])

        x=coord[0]
        y=coord[1]
        
        
        #print("\n",self.matrix_backtracking)
        #print("\n",self.matrix_scores)
        self.recursion_back(x, y, (x,y)) # i call the recursion back
        #print("\n\n Dictionary:\n",self.dict_sequences)

        #self.alignment_scores() # i call teh function that will cfreate the alignments


    def recursion_back(self,i, j, score): # i = rows, j = columns

        if type(self.matrix_backtracking[i,j]) == str:
            list_movements=self.matrix_backtracking[i,j].split("\t")
        else:
            list_movements=self.matrix_backtracking[i,j].decode("utf-8").split("\t") # i transform in a list the string that is contained in the backtracking matrix with the symbols for the provenience of the score (D for diagonal, L for left, U for up)
        
        #print(list_movements)
        #print(self.list_alig)
        if len(list_movements) == 1:
            #print("lunghezza list pari ad 1")
    
            if list_movements[0]=="D":
                self.list_alig[self.temp]+="D"
                #print("Align:",self.list_alig)
                self.recursion_back(i-1, j-1,score)

            elif list_movements[0]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[0]=="U":   
                self.list_alig[self.temp]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)
            elif list_movements[0]=='':
                if score not in self.dict_sequences:
                    self.dict_sequences[score]=[self.list_alig[self.temp]]
                else:
                    self.dict_sequences[score].append(self.list_alig[self.temp])
                self.list_alig.remove(self.list_alig[self.temp])
                self.temp=len(self.list_alig)-1
            
            #print(self.list_alig)
            #print(self.temp)
        
        if len(list_movements) == 2:
            #print(list_movements[0], list_movements[1])
            # con temp

            self.list_alig.append(self.list_alig[self.temp])
            
            #print(self.list_alig)
            
            #print(list_movements)
            if list_movements[0]=="D":
                    
                self.list_alig[self.temp]+="D"
                #print(self.list_alig)
                    
                self.recursion_back(i-1, j-1,score)
            elif list_movements[0]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[0]=="U":
                    
                self.list_alig[self.temp-1]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)
            
                
            

            #self.listina.append(self.listina[self.temp-1][:-2])
            #print(self.list_alig)

            if list_movements[1]=="D":
                    
                self.list_alig[self.temp]+="D"
                #print(self.list_alig)
                    
                self.recursion_back(i-1, j-1,score)
            elif list_movements[1]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[1]=="U":
                    
                self.list_alig[self.temp]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)
            
            #print(self.list_alig)
            #print(self.temp)

        if len(list_movements) == 3:
            # con temp

            self.list_alig.append(self.list_alig[self.temp])
            self.list_alig.append(self.list_alig[self.temp])
            #print(self.list_alig)

            if list_movements[0]=="D":
                    
                self.list_alig[self.temp]+="D"
                #print(self.list_alig)
                    
                self.recursion_back(i-1, j-1,score)
            elif list_movements[0]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[0]=="U":
                    
                self.list_alig[self.temp]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)

            #print(self.list_alig)

            if list_movements[1]=="D":
                    
                self.list_alig[self.temp]+="D"
                #print(self.list_alig)
                    
                self.recursion_back(i-1, j-1,score)
            elif list_movements[1]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[1]=="U":
                    
                self.list_alig[self.temp]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)


            #print(self.list_alig)

            if list_movements[2]=="D":
                    
                self.list_alig[self.temp]+="D"
                #print(self.list_alig)
                    
                self.recursion_back(i-1, j-1,score)
            elif list_movements[2]=="L":
                self.list_alig[self.temp]+="L"
                #print(self.list_alig)
                self.recursion_back(i, j-1,score)
            elif list_movements[2]=="U":
                    
                self.list_alig[self.temp]+="U"
                #print(self.list_alig)
                self.recursion_back(i-1, j,score)
            
            #print(self.dict_sequences)
            #print(self.temp)
 

    def alignment_scores(self): ## NB the sequences are reversed! I need to reverse them!
        alignments_with_score={}
        seq1=self.getSeq1()
        seq2=self.getSeq2()
        for key in self.dict_sequences: # i iterate over the keys of the dict -> meaning all the starting positions of all the paths i have

            score=self.matrix_scores[key]
            for path in self.dict_sequences[key]: # iterate over all the possible paths associated to the initial position
                i=key[0]
                j=key[1]
                #print(path)
                sequences=["","",""]
                #print(sequences)
                for direction in path:
                    if direction=="D":
                        sequences[0]+=seq1[i-1]
                        sequences[1]+=seq2[j-1]
                        if seq1[i-1]==seq2[j-1]:
                            sequences[2]+="M"
                        else:
                            sequences[2]+="S"
                        i=i-1
                        j=j-1
                    if direction=="U":
                        sequences[0]+=seq1[i-1]
                        sequences[1]+="_"
                        sequences[2]+="G"
                        i=i-1
                    if direction=="L":
                        sequences[0]+="_"
                        sequences[1]+=seq2[j-1]
                        sequences[2]+="G"
                        j=j-1
                if score not in alignments_with_score:
                    alignments_with_score[score]=[]
                alig_reverse=["","",""]
                alig_reverse[0]=sequences[0][::-1]
                alig_reverse[1]=sequences[1][::-1]
                alig_reverse[2]=sequences[2][::-1]
                alignments_with_score[score].append(alig_reverse)
            
            #print(alignments_with_score)
            
        self.print_results(alignments_with_score)        


    def print_results(self, dictionary_score):

        # lets print the best result with the highest score:
        max_now=0
        
        

        for key in dictionary_score:
            print(key)
            max_now=max(max_now, key) # i find the max score looking at the position of teh keys in the dictionary
            
            
        print("ALIGNMENT RESULTS:\t LEGEND SYMBOLS: *: match; \t |: mismatch; \n")
        print("SEQUENCE 1: \t", self.getSeq1())
        print("SEQUENCE 2: \t", self.getSeq2())
        print("\n")
        print("Score max:\t", max_now)
        print("\n")

        
        for value in dictionary_score[max_now]:
            matches=0
            mismatches=0
            gaps=0
            print(value[0])
            stringhetta=""
            for element in value[2]: # i iterate over the scores:
                if element=="M": # match
                    stringhetta=stringhetta+"*"
                    matches=matches+1
                elif element=="S": # mismatch
                    stringhetta=stringhetta+"|"
                    mismatches=mismatches+1
                else: # i have a gap
                    stringhetta=stringhetta+" "
                    gaps=gaps+1
            print(stringhetta)
            print(value[1])
            print("\nAdditional info:")
            print("Length of the sub-alignment: ", len(value[2]), "\t Total mismatches: ",mismatches, "\t Total matches:", matches, "\t Total gaps: ", gaps)

            print("\n"+50*"-"+"\n")


#### Example 1: The modified version of the script must return all the possible local alignments, ordered by decreasing score, with length > 5, score > 4 and gaps > 0. Gaps are included in the calculation of alignment length. Overlaps in trace-back paths are allowed. 

### The modified version of the program must return all the possible local alignments (can be more than one) with the minimum number of gaps (can be 0) among the ones that have a length of at least 7. Gaps are included in the calculation of alignment length. Overlaps in trace-back paths are allowed.



if __name__ == '__main__': # for the code to be called from the terminal

    # creation of the Parser that will be needed to write user-friendly command-line interfaces -> call the command SmithWaterman 
    # use -h to get the helper, not needed to explicit it, there is already in default

    parser = ap.ArgumentParser(prog='Smith Waterman algorithm',  description='Smith Waterman command to perform a local alignment', epilog='This was the help!  The output will be written in a file.txt Thanks!')
    parser.add_argument('seq1', type= str, help='First of the two string to align')
    parser.add_argument('seq2', type= str, help='Second string to align')
    parser.add_argument('-m', '--match', type=float, default=3, help='positive number corresponding to the score for a match in the sequence alignment. [default: 3]')
    parser.add_argument('-mm', '--mismatch', type=float, default=-2, help='negative number corresponding to the score for a mismatch in the sequence alignment. [default: -2]')
    parser.add_argument('-g', '--gap', type=float, default=-1, help='negative number corresponding to the score for a gap in the sequence alignment. [default: -1]')
    args = parser.parse_args()

    

    sw=SmithWaterman(args.seq1, args.seq2, args.match, args.mismatch, args.gap) # calling the class SmithWaterman and giving in input the values defined by the user (or in default if not given, for the optional paramerers)
    sw.M_Filling()
    sw.greedy()
    sw.alignment_scores()
