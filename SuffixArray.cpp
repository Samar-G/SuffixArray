#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include<chrono>

using namespace std;

//define the director where the file we need to read lies
#define DIR "E:\\Bioninformatics\\year3,sem2\\Sequence Analysis\\Assignments\\"

//used to read the first 10k characters
char* readFile()
{
    //open the file and declaring the we need to read from it
    FILE* file=fopen(DIR "genome.fasta", "r");
    char* buf=new char[10001];

    buf[0]=0;
    //To ignore the id line(first line in each genome)
    fscanf(file, "%[^\n\r] ", buf);

    //to make sure it doesn't read over something already there
    buf[0]=0;
    //to read the first 10k characters
    fscanf(file, "%10000[^\r] ", buf);
    //to print the characters(string) to screen
    printf("%s\n", buf);
    //return the character array so I can use it later in the main
    return buf;
}

//used to read the whole genomes
char* readFileWhole()
{
    FILE* file=fopen(DIR "genome.fasta", "r");
    char* buf=new char[3000001];

    buf[0]=0;
    //To ignore the id line(first line in each genome)
    fscanf(file, "%[^\n\r] ", buf);

    //to make sure it doesn't read over something already there
    buf[0]=0;
    //to read the first 10k characters
    fscanf(file, "%2913956[^\r] ", buf);
    //ro print the characters(string) to screen
    printf("%s\n", buf);
    //return the character array so I can use it later in the main
    return buf;
}

//global of the sting that I will be sorting to make it easier for access across the code
char* seq;
//intializes global size to be used later due to some struggles faced for sorting and will be assigned a number in the main through each test case
int n;

/* in this function we compare each suffix by refereing to the main char array(string)
   checking if the i and j are within the size as once it is not we will break
   because one will just end up comparing with nothing if the other is out of size.
   */
int compareFirst(int i, int j)
{
    do
    {
        if (seq[i] > seq[j])
        {
            return 0;
        }
        else if (seq[i] < seq[j])
        {
            return 1;
        }
        i++;
        j++;
    } while (i<n && j<n);
    if(j == n)
    {
       return 0;
    }
    return 1;
}


/* through this function we build the suffix array using the naive method O(n^2 logn)
   first we initialize an array of size n and inside is the ordered from 0 to size-1 and then it will be sorted
   using the built in and a comparator we created above*/
int* buildNaive()
{
    int* suffix= new int[n];

    //initializing from 0 to size-1
    for(int i=0; i<n; i++)
    {
        suffix[i]=i;
    }
    //sorting accordingly to the comparator
    sort(suffix,suffix+n,compareFirst);

    return suffix;
}


//writing to the file, takes the arra or the pointer to array that has the thing we need to write
int* writeFile(int* buf)
{
    FILE* file=fopen("out.txt", "w");
    //loops over to write what is inside
    for(int i=0; i<n; i++)
    {
        fprintf(file,"%d",buf[i]);
    }
    fclose(file);
    return buf;
}

int main()
{
    seq=readFile(); //reads first 10k chars
    n=100001;
    auto beg = std::chrono::steady_clock::now(); //checks the time before calling the nave function
    int *tenkChars=buildNaive(); //calls th enaive function
    writeFile(tenkChars);   //writes the suffix array to a file
    auto fin = std::chrono::steady_clock::now(); //check the time after done with function
    std::chrono::duration<double> seconds = fin-beg; //calculate time in seconds
    cout<<seconds.count() << "s\n"; //outputs total time after run in second

    //read the whole genome
    seq=readFileWhole();
    n=3000000;
    int *whole=buildPrefix();
    writeFile(whole);


/*SOME TEST CASES*/

    seq= {"ACTGAGTGTAGTA"};
    n=13;
    int *se13=buildNaive();

    for (int i = 0; i < n; i++)
    {
        cout << se13[i] << " ";
    }

    seq= {"AATTAGGGTAGTG"};
    n=12;
    int *se12=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se12[i] << " ";
    }

    seq= {"TATTTTCCCTTTTATAAGTTT"};
    n=21;
    int *se11=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se11[i] << " ";
    }

    seq= {"GTGATGGTCTTTGGTTTTCTATA"};
    n=23;
    int *se10=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se10[i] << " ";
    }

    seq= {"GTTAAAAAATTAG"};
    n=13;
    int *se=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se[i] << " ";
    }

    seq= {"AAAAAAAAAAA"};
    n=11;
    int *se1=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se1[i] << " ";
    }

    seq= {"AAACAGGGTCCCT"};
    n=13;
    int *se2=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se2[i] << " ";
    }

    seq= {"TATTTGAAGTTACAT"};
    n=15;
    int *se3=buildNaive();

    for (int i = 0; i < n; i++)
    {
        cout << se3[i] << " ";
    }

    seq= {"CGTTATAGCAACTA"};
    n=14;
    int *se4=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se4[i] << " ";
    }

    seq= {"AACATCACCATGAGTT"};
    n=16;
    int *se5=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se5[i] << " ";
    }

    seq= {"GTGATGGTCTTTTTCTATA"};
    n=19;
    int *se6=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se6[i] << " ";
    }

    seq= {"GCAAGATGAGCAA"};
    n=12;
    int *se7=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se7[i] << " ";
    }

    seq= {"AAAGTCAAAAA"};
    n=11;
    int *se8=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se8[i] << " ";
    }

    seq= {"ACCTATATAATTAT"};
    n=14;
    int *se9=buildNaive();
    for (int i = 0; i < n; i++)
    {
        cout << se9[i] << " ";
    }

    return 0;
}
