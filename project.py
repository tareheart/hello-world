# species named in the data file finches.csv 
finch1 = 'G.conirostris_Espanola'
finch2 = 'G.conirostris_Genovesa'
finch3 = 'G.difficilis'
finch4 = 'G.magnirostris'
outgroup = 'L.noctis'
            

def read_input(fileName):
    ''' The goal of this function is to take one argument which will be a string name of a csv file, then to create a list combining the A and B allele of an individual that shares the same id. It will grab necessary items from row 1 as well as row 2, ultimately appending the combined information as a list of a list for all rows in the file. '''
    import csv
    fileref = open(fileName, "r")  # opens the file in reading mode
    data_reader = csv.reader(fileref)
    inner_list = []  # sets up inner list and final list in order to simplify future coding
    final_list = []

    for row in data_reader:  # iterates over each row in the finches file
        if row[2] == "A": # identifies the 1st row (A allele) of the individual and adds necessary items
            inner_list.extend([row[0], row[1], row[3], None, float(row[4]), None, row[5]])
        elif row[2] == "B":  # identifies the 2nd row of the same individual
            if inner_list:  # if statement to make sure that "inner_list" already exists
                inner_list[3] = row[3]  # replaces the None in the existing list with necessary item
                inner_list[5] = float(row[4])  # does the same but replaces it as a float value
                final_list.append(inner_list)  # adds the combined inner list to final_list
                inner_list = []  # "empties out" inner_list to add new/different items
    return final_list  # after for loop is finished iterating, returns final_list
    


def allele_dist(gene1, gene2):
    ''' This function takes in two arguments which are both strings. The goal is to calculate and return the hamming distance which is essentially the number of positions at which the two strings are different. '''
    ham_dist = 0  # assigns 0 to a variable called ham_dist

    # assumes that the two strings are of the same distance, hence the len(gene1)
    for allele in range(0, len(gene1)):  # for loop iterates in the range of 0 to length of gene
        if gene1[allele] != gene2[allele]:  # adds 1 to ham_dist for every symbol that is different
            ham_dist += 1

    return ham_dist  # returns the hamming distance
        


def gene_dist(finch1, finch2):
    ''' For this function, two argument which are both lists are taken. The goal is to compute the overall gene distance using the allele_dist function and averaging the result out. '''
    dist1 = allele_dist(finch1[2], finch2[2])  # calling to function allele_dist with given 4 sets of data
    dist2 = allele_dist(finch1[2], finch2[3])  # repeated four times and each assigned to a variable
    dist3 = allele_dist(finch1[3], finch2[2])
    dist4 = allele_dist(finch1[3], finch2[3])

    avg_dist = (dist1 + dist2 + dist3 + dist4) / 4  # ultimately added and divided by 4 to find the average

    return avg_dist  # the avg_dist is returned



def beak_dist(finch1, finch2):
    ''' This function takes in 2 arguments that are both lists. The goal is to calculate the difference between beak distance in the two lists, find the average and to ultimately return it. '''
    dist1 = abs(finch1[4] - finch2[4])  # uses the four sets of data but locates for the beak distance
    dist2 = abs(finch1[4] - finch2[5])  # uses the abs() for any negative values that may occur
    dist3 = abs(finch1[5] - finch2[4])  # e.g |3-6| = 3 rather than -3
    dist4 = abs(finch1[5] - finch2[5])  # assigns each calculation to 4 variables

    avg_dist = (dist1 + dist2 + dist3 + dist4) / 4  # each variable is added and divided by 4 

    return avg_dist  # the avg_dist is returned as a float 



def outgroup_distance(finches, speciesName, outgroupName):
    ''' This function will take in three arguments, where first is the finches list returned by the read_input(function), second is the scientific name of Darwin's species, and third is the scientific name for the outgroup species. The goal of this function is to compute gene distance and beak distance from the outgroup to each Darwin's species. Then the values are added to an existing list and both lists are eventually returned as a tuple. '''
    geneDist = []  # creates an empty list for gene distance
    beakDist = []  # creates an empty list for beak distance
    non_out = []  # creates an empty list for Darwin species 

    for sublist in finches:  # for loop to separate the list for outgroup and Darwin's species
        if speciesName in sublist:  # looks for Darwin's species name and adds the list if it is found
            non_out.append(sublist)
        elif outgroupName in sublist:  # looks for outgroup species and assigns it as outlist if found
            outlist = sublist

    for sublist in non_out:  # for loop iterates over each list inside Darwin species' list of lists
        geneDist.append(gene_dist(sublist, outlist))  # calls gene_dist with sublist and outlist
        beakDist.append(beak_dist(sublist, outlist))  # adds to the appropriate list depending on beak or gene
    return geneDist, beakDist  # returns both lists as a tuple



''' Add code below to call/use/verify above functions and plot results '''
import matplotlib.pyplot as plt

finches = read_input("finches.csv")
specie1 = outgroup_distance(finches, finch1, outgroup)  # call outgroup_distance function for 4 groups
specie2 = outgroup_distance(finches, finch2, outgroup)  # these will return a tuple of geneDist, beakDist
specie3 = outgroup_distance(finches, finch3, outgroup)  # contain multiple values for each individual
specie4 = outgroup_distance(finches, finch4, outgroup)

x1 = specie1[0]  # assign appropriate list. [0] will be assigned to list of geneDist (with multiple values)
y1 = specie1[1]  # [1] will be assigned to list of beakDist (with multiple values)
x2 = specie2[0]
y2 = specie2[1]
x3 = specie3[0]
y3 = specie3[1]
x4 = specie4[0]
y4 = specie4[1]

plt.scatter(x1, y1, c='b', label=finch1)  # 4 different plot statements that create a scatter plot
plt.scatter(x2, y2, c='r', label=finch2)  # each plot statement uses a different color
plt.scatter(x3, y3, c='g', label=finch3)  # label created for legend
plt.scatter(x4, y4, c='y', label=finch4)

plt.xlabel('gene distance to the outgroup')  # legend, title, and x/y labels are added
plt.ylabel('beak distance to the outgroup')
plt.title("Darwin's Finches vs Outgroup Finch Data")
plt.legend()

plt.savefig("Plot1.png")
