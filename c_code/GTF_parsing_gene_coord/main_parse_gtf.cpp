#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
//#include <algorithm>
#include <string>
#include <iostream>

#include "functions_parse_gtf.hpp"

static struct option long_opts[]={
    /* ... */
    {"gtf_file", required_argument, NULL, 'g'}, //use this for annotation input file
    {"out_file", required_argument, NULL, 'o'}, //use this for output file
    {NULL, 0, NULL, 0}};
static char short_option[]="g:o:h"; //use this for parsing // this could have  *short_option

int main(int argc, char *argv[])
{
    
    char *gtf_file = NULL;    // arg for GTF file .
    char *output_file = NULL; // arg for output  file .

    assign_variables(argc, argv, short_option, long_opts, &gtf_file, &output_file);
    
    if (gtf_file == NULL)
    {
        printf("We do not have GTF file. Exiting\n");
        exit(0);
    }

    if (output_file == NULL)
    {
        printf("We do not have output file. Exiting\n");
        exit(0);
    }

    if (argc != 5)
    {
        printf("we do not have sufficient arguments\n");
        print_usage();
        exit(0);
    }

    FILE *fp_gtf = fopen(gtf_file, "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf     //FILE *fp_gtf = fopen("test.gtf", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf

    if (fp_gtf == NULL)
        exit(EXIT_FAILURE);

    char *line = NULL;
    size_t len = 0;

    const char *delimiters = "\t";                                           //use this for spliting lines
    std::map<std::string, std::pair<std::string, std::string>> map_gene_chr; // APOE(key) APOE CHR19 <value> (pair)

    while ((getline(&line, &len, fp_gtf)) != -1)
    {
        if (line[0] != '#')
        {
            gene gene_temp_info = gene(); //make a new struct of gene
            char *token = NULL;
            remove_trailingspaces(line); //remove trailing new line character

            token = strtok(line, delimiters);
            int count_split = 1; //use this for counting split. as columns are fixed
            int found_gene = 0;  //if column 3rd is gene             //char gene_info_print[5000]; //use this to print CHR1 start end ENS_ID GENE_NAME origin and biotype

            while (token != NULL)
            {
                if (count_split == 1)
                {
                    if (check_chromosome(token) == 0)
                    {
                        break;
                    }
                    else
                    {
                        gene_temp_info.chromosome = token; //mystrcat(gene_info_print, "chr");  mystrcat(gene_info_print, token);
                    }
                }
                if (count_split == 3 && (strcmp(token, "gene") == 0))
                {
                    found_gene = 1;
                }
                if (count_split == 3 && (strcmp(token, "gene") != 0))
                {
                    break;
                }
                if (count_split == 4 && found_gene == 1)
                {
                    /**
                     * join with start position
                    */
                    gene_temp_info.start_position = token; //mystrcat(gene_info_print, "\t");  //mystrcat(gene_info_print, token); //printf("%s\n", token);
                }

                if (count_split == 5 && found_gene == 1)
                {
                    /**
                     * join with end position
                    */
                    gene_temp_info.end_position = token; //mystrcat(gene_info_print, "\t"); //mystrcat(gene_info_print, token);
                }

                if (count_split == 9 && found_gene == 1)
                {
                    char *gene_info = NULL;
                    gene_info = strtok(token, ";"); // gene_id "ENSG00000223972";gene_name "DDX11L1";gene_source "ensembl_havana";gene_biotype "pseudogene";

                    while (gene_info)
                    {
                        if (gene_info[5] == 'i') // gene_id
                        {
                            remove_double_quotes(&gene_info); //gene_id "ENSG00000223972";as ENSG00000223972

                            gene_temp_info.gene_id = gene_info; // mystrcat(gene_info_print, "\t");  // mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 'n') //gene_name
                        {
                            remove_double_quotes(&gene_info);     //make gene_name "DDX11L1"  as DDX11L1
                            gene_temp_info.gene_name = gene_info; // printf("the name of genename ided is %s\n", gene_info);//  mystrcat(gene_info_print, "\t"); //mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 's') //  gene_
                        {
                            remove_double_quotes(&gene_info);       //make gene_source "ensembl_havana" as  ensembl_havana //printf("the name of source ided is %s\n", gene_info);
                            gene_temp_info.GENE_SOURCE = gene_info; //  mystrcat(gene_info_print, "\t"); //mystrcat(gene_info_print, gene_info);
                        }

                        if (gene_info[6] == 'b') //gene_biotype
                        {
                            remove_double_quotes(&gene_info);        //make gene_biotype "pseudogene" as  pseudogene
                            gene_temp_info.GENE_BIOTYPE = gene_info; // printf("the name of biotype ided is %s\n", gene_info); //mystrcat(gene_info_print, "\t"); //mystrcat(gene_info_print, gene_info);
                        }

                        gene_info = strtok(NULL, ";");
                    }

                    // check if gene exists. If it exists then add _dup to its name

                    if (map_gene_chr.count(gene_temp_info.gene_name) > 0)
                    {
                        std::cout << gene_temp_info.gene_name << " is an element of map\n";

                        if ((map_gene_chr[gene_temp_info.gene_name].second).compare(gene_temp_info.chromosome) == 0)
                        {
                            //check if same CHR. Then only add _dup. Same gene name on diff CHR aren't concern at the moment
                            gene_temp_info.gene_name = gene_temp_info.gene_name + "_dup";
                            map_gene_chr.insert({gene_temp_info.gene_name, {gene_temp_info.gene_name, gene_temp_info.chromosome}});
                        }
                        else
                        {
                            std::cout << gene_temp_info.gene_name << " same gene name but Diff CHR\n";
                        }
                    }
                    //if ends if gene already exists

                    else
                    {
                        //gene name is new and has not been added to the map                         //https://stackoverflow.com/a/59484450/2740831
                        map_gene_chr.insert({gene_temp_info.gene_name, {gene_temp_info.gene_name, gene_temp_info.chromosome}});
                    }

                    std::cout << gene_temp_info.chromosome << '\t' << gene_temp_info.gene_id << '\t' << gene_temp_info.gene_name
                              << '\t' << gene_temp_info.start_position << '\t'
                              << gene_temp_info.end_position << '\t'
                              << gene_temp_info.GENE_SOURCE << '\t'
                              << gene_temp_info.GENE_BIOTYPE << '\t' << '\n';

                    break; //move out of the while loop of line parsing.
                }
                ///parsing of column with gene info ends

                token = strtok(NULL, delimiters);
                count_split++;
            }
            //while loop ends for token split of delimiters
        }
        //if line doesn't start with an #

    } //while loop of file ends

    //while loop ends of file reading
    if (line)
    {
        free(line);
    }

    fclose(fp_gtf); //close file

    std::map<std::string, std::pair<std::string, std::string>>::iterator it_map_gene_pair;

    for (it_map_gene_pair = map_gene_chr.begin(); it_map_gene_pair != map_gene_chr.end(); it_map_gene_pair++)
    {
        std::cout << it_map_gene_pair->first << '\t' << (it_map_gene_pair->second).first << "\t"
                  << (it_map_gene_pair->second).second << '\n'; //
    }
    return 0;
}
//main function ends

// g++ c14 supported https://stackoverflow.com/a/34681870/2740831
/**
     * take in GTF file. Print on the terminal CHR start, end, ensembl id, gene name and gene-biotype
     * //  gcc -Wpedantic -Wextra -Wall main_parse_gtf.c -o gtf_parse
     * //  make clean ; make ; ./gtf_parse  > genes
     * //  cat test.gtf | grep gene | grep -v exon | grep -v transcript
     * //  printf "CHR    START    END    ENSEMBL_IDGENE_NAME    GENE_SOURCE    GENE_BIOTYPE\n" >genes_CHR37
    */

//FILE *fp_gtf = fopen("/mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.gtf", "r"); // /mnt/mfs/ctcn/resources/GRCh38/v1/Homo_sapiens.GRCh38.93.filtered.gtf