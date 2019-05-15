import pandas as pd
import urllib as url
import csv
import os
import xml.etree.cElementTree as et

base_path = os.path.dirname(os.path.realpath(__file__))
XMLcurrProteinID = ""

targetType = []

BioAssayDF = pd.DataFrame()
nanTargeturlCounter = 0
NOTnanCounter = 0

nanPubChemCIDCounter = 0
NOTnanPubChemCIDCounter = 0

yesUniprotFASTA = 0
noUniprotFASTA = 0

outputFile = open("./outputs.txt", 'w')

drugDF= pd.DataFrame()
drugDF = pd.read_csv('drugDatasetDrugBank.txt', delimiter='\t', low_memory=False, header= 0,
                                      names=["drugbankId", "drugName", "drugSMILES", "drugInChI", "drugPubChemCompound" ])

fromList = ["ACC","ID","UPARC","NF90","NF100","GENENAME","EMBL_ID","EMBL","P_ENTREZGENEID","P_GI",
            "PIR","REFSEQ_NT_ID","P_REFSEQ_AC","UNIGENE_ID","PDB_ID","DISPROT_ID","BIOGRID_ID",
            "DIP_ID","MINT_ID","STRING_ID","CHEMBL_ID","DRUGBANK_ID","GUIDETOPHARMACOLOGY_ID",
            "SWISSLIPIDS_ID","ALLERGOME_ID","ESTHER_ID","MEROPS_ID","MYCOCLAP_ID","PEROXIBASE_ID",
            "REBASE_ID","TCDB_ID","BIOMUTA_ID","DMDM_ID","WORLD_2DPAGE_ID","DNASU_ID","ENSEMBL_ID",
            "ENSEMBL_PRO_ID","ENSEMBL_TRS_ID","ENSEMBLGENOME_ID","ENSEMBLGENOME_PRO_ID",
            "ENSEMBLGENOME_TRS_ID","GENEDB_ID","P_ENTREZGENEID","KEGG_ID","PATRIC_ID","UCSC_ID",
            "VECTORBASE_ID","WBPARASITE_ID","ARACHNOSERVER_ID","ARAPORT_ID","CCDS_ID","CGD",
            "CONOSERVER_ID","DICTYBASE_ID","ECHOBASE_ID","ECOGENE_ID","EUHCVDB_ID","EUPATHDB_ID",
            "FLYBASE_ID","GENECARDS_ID","GENEREVIEWS_ID","H_INVDB_ID","HGNC_ID","HPA_ID",
            "LEGIOLIST_ID","LEPROMA_ID","MAIZEGDB_ID","MGI_ID","MIM_ID","NEXTPROT_ID","ORPHANET_ID",
            "PHARMGKB_ID","POMBASE_ID","PSEUDOCAP_ID","RGD_ID","SGD_ID","TUBERCULIST_ID","WORMBASE_ID",
            "WORMBASE_PRO_ID","WORMBASE_TRS_ID","XENBASE_ID","ZFIN_ID","EGGNOG_ID","GENETREE_ID",
            "HOGENOM_ID","HOVERGEN_ID","KO_ID","OMA_ID","ORTHODB_ID","TREEFAM_ID","BIOCYC_ID",
            "REACTOME_ID","UNIPATHWAY_ID","CLEANEX_ID","COLLECTF_ID","CHITARS_ID","GENEWIKI_ID",
            "GENOMERNAI_ID"]

for index, row in drugDF.iterrows():

    currPubChemCID = row["drugPubChemCompound"]

    if(str(currPubChemCID)!="nan"):

        #print("** currPubChemCID: " + str(currPubChemCID))

        NOTnanPubChemCIDCounter = NOTnanPubChemCIDCounter + 1

        urlLink = "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=jsonp&query=%5B%7B%22download%22:%5B%22activity%22,%22acvalue%22,%22acname%22,%22targetname%22,%22targeturl%22,%22aidname%22,%22aid%22,%22sid%22,%22cid%22%5D,%22collection%22:%22bioactivity%22,%22where%22:%7B%22ands%22:%5B%7B%22cid%22:%22" + str(currPubChemCID) + "%22%7D%5D%7D,%22order%22:%5B%22acvalue,asc%22%5D,%22start%22:1,%22limit%22:1000000,%22nullatbottom%22:1%7D,%7B%22histogram%22:%22activity%22,%22bincount%22:10000%7D%5D"

        #BioAssayDF = url.urlopen(urlLink).read()
        BioAssayDF = pd.read_table(urlLink,sep=",", index_col= 0)

        # extract the target type
        # for each row in BioAssayDF
        for indexBioAssay, rowBioAssay in BioAssayDF.iterrows():

            currTargeturl = rowBioAssay["targeturl"]

            #print(currTargeturl)

            # if NULL the type is float (I don't know why)
            if (str(currTargeturl)=="nan"):
                nanTargeturlCounter = nanTargeturlCounter + 1

                currTargetName = rowBioAssay["targetname"]
                if(str(currTargetName)!="nan"):
                    print("###" + currTargetName)
                    outputFile.write(currTargetName)

            else:

                FASTAFlag = False

                NOTnanCounter = NOTnanCounter + 1
                if(currTargeturl[0:8]=="/target/"):
                    nextBackSlash = currTargeturl[9:len(currTargeturl)].find('/')

                    # target type
                    currTargetType = currTargeturl[8:(nextBackSlash + 9)]
                    # we need only protein type
                    if (currTargetType == "protein"):

                        currProteinID = currTargeturl[(nextBackSlash + 10): len(currTargeturl)]

                        print(str(currPubChemCID) + " *** " + currProteinID )
                        outputFile.write(str(currPubChemCID) + " *** " + str(currProteinID))

                        # get the another proteinID from the /target/protein/???.xml
                        htmlProteinURL = "https://pubchem.ncbi.nlm.nih.gov/target/protein/" + str(currProteinID)

                        # Error after: 9878 *** ABY84639 *** ABY84639
                        try:
                            htmlProteinFileTxt = url.urlopen(htmlProteinURL).read()  # error !!!
                        except IOError:
                            print("IOError")
                            pass

                        #keyString = "https://pubchem.ncbi.nlm.nih.gov/target/protein"
                        keyString = "data-pubchem-id="
                        startIndex = htmlProteinFileTxt.find(keyString)
                        endIndex = startIndex + len(keyString)
                        proIDend = htmlProteinFileTxt.find('"', endIndex+1)
                        currAltProteinID = htmlProteinFileTxt[endIndex + 1:proIDend]

                        #print(str(currPubChemCID)+ " *** " + currProteinID + " *** " + currAltProteinID)

                        if (currAltProteinID==""):
                            print("Err2 : " + currProteinID + " *** " + currAltProteinID)
                            outputFile.write("Err2 : " + currProteinID + " *** " + currAltProteinID)

                        uniprotURL = "https://www.uniprot.org/uniprot/" + str(currAltProteinID) + ".fasta"
                        try:
                            uniprotFASTA = url.urlopen(uniprotURL).read()
                        except IOError:
                            print("IOError")
                            pass

                        if((len(uniprotFASTA)>0 and uniprotFASTA[0]=="<") or (len(uniprotFASTA)==0)):

                            # ----------------------------------
                            # try to get the uniprotKB ID using

                            currAltUniPortFlag = False
                            base = 'http://www.uniprot.org'
                            tool = 'mapping'

                            for fromType in fromList:
                                params = {'from': fromType,
                                          'to': "",
                                          'format': 'tab',
                                          'query': currAltProteinID, # may currProteinID ?????
                                          }
                                data = url.urlencode(params)
                                mapURL = base + '/' + tool + '?' + data

                                #error: 6106.0 *** 1Y7V_A
                                response = url.urlopen(mapURL)

                                try:
                                    output = response.read() #provides tab-delimited output of the mapping
                                except IOError:
                                    print("IOError")
                                    pass

                                if (len(output) > 8):
                                    output = output.replace('\n', '\t')
                                    outputs = output.split('\t')
                                    currAltUniPort = outputs[3]
                                    currAltUniPortFlag = True
                                    print("$$ " + fromType + " $$" + currAltUniPort)
                                    outputFile.write("$$ " + currAltUniPort)

                            # ----------------------------------

                            if(currAltUniPortFlag==True):

                                print("FASTA")
                                outputFile.write("FASTA")
                                yesUniprotFASTA = yesUniprotFASTA + 1
                                uniprotURL = "https://www.uniprot.org/uniprot/" + str(currAltUniPort) + ".fasta"
                                uniprotFASTA = url.urlopen(uniprotURL).read()
                                FASTAFlag = True

                            else:

                                print("not FASTA")
                                outputFile.write("not FASTA")
                                noUniprotFASTA = noUniprotFASTA + 1

                                xmlProteinURL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/protein/" + currAltProteinID + \
                                            "/XML/?response_type=save&response_basename=PROTACXN_" + currAltProteinID
                                xmlProteinFileName = "./xmlProtein/" + currProteinID + ".xml"

                                xmlProteinFile = open(xmlProteinFileName, 'w')
                                xmlProteinFile.write(url.urlopen(xmlProteinURL).read())
                                xmlProteinFile.close()

                                xmlProtein = os.path.join(base_path, "./xmlProtein/" + currProteinID + ".xml")
                                tree = et.parse(xmlProtein)
                                root = tree.getroot()

                                for child in root:  # root is <Record> and the child
                                    # print("1. " + child.tag)
                                    for nextChild in child:

                                        if ("Information" in nextChild.tag):
                                            for FASTAnextChild in nextChild:
                                                if ("StringValue" in FASTAnextChild.tag):
                                                    pubchemFASTA = FASTAnextChild.text
                                                    #print(pubchemFASTA)
                                                    #print("+++++++++++++++++++++++++++++++++++++++++")

                                                    if(pubchemFASTA!=""):
                                                        FASTAFlag = True

                        elif(len(uniprotFASTA)>0):

                            print("yes FASTA")
                            outputFile.write("yes FASTA")

                            yesUniprotFASTA = yesUniprotFASTA + 1
                            #print("&& uniprotFASTA &&")
                            #print(uniprotFASTA)
                            #print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

                            if (uniprotFASTA != ""):
                                FASTAFlag = True

                        else:
                            print("No FASTA " + str(currProteinID) + " !!!!!!!! " )
                            outputFile.write("No FASTA " + str(currProteinID) + " !!!!!!!! " )

                        if (FASTAFlag==False):
                            print("Err: " + currProteinID + " *** " + str(currPubChemCID))
                            outputFile.write("Err: " + currProteinID + " *** " + str(currPubChemCID))

        #print("-------------------------------")

    else:

        nanPubChemCIDCounter = nanPubChemCIDCounter + 1

#proteinsDataset.to_csv("proteinsDataset.txt", sep="\t", index=False)

print(nanPubChemCIDCounter)
print(NOTnanPubChemCIDCounter)
print(nanTargeturlCounter)
print(NOTnanCounter)
print("----------------------")
print(yesUniprotFASTA)
print(noUniprotFASTA)

outputFile.close()