#define _CRT_SECURE_NO_WARNINGS // a Microsoft thing about strcpy security issues

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <stdlib.h>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <vector>
#include <utility>
#include <dirent.h>
#include <omp.h>
#include <string.h>
#include "./fasta_reader/fasta.h"

using namespace std;

// graceful exit of sorts
void panic(const char *msg)
{
    cerr << msg << endl;
    exit(0);
}

// approximative timing of various steps
struct tTimer
{
    time_t begin;
    tTimer() { begin = time(NULL); }
    void print(const char *msg)
    {
        time_t now = time(NULL);
        cerr << msg << " in " << now - begin << "s\n";
        begin = now;
    }
};

// load input pattern
void read_sequence(const char *filename, long int &len, int &num_contigs)
{
    //Fasta reader part
    FASTAFILE *ffp = OpenFASTA(filename);
    char *seq;
    char *name;
    int L;
    int contigs = 0;

    long int total_length = 0;
    while (ReadFASTA(ffp, &seq, &name, &L))
    {
        total_length += L;
        contigs++;

        free(seq);
        free(name);
    }

    len = total_length;
    num_contigs = contigs;

    CloseFASTA(ffp);
}

// read all FASTA files in directory
vector<string> read_directory(char *pathname)
{
    DIR *dir;
    struct dirent *ent;
    string path_name(pathname);

    vector<string> path_to_files;

    if ((dir = opendir(pathname)) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL)
        {
            if (!strcmp(ent->d_name, ".") or !strcmp(ent->d_name, ".."))
                continue;
            path_to_files.push_back(path_name + ent->d_name);
        }

        closedir(dir);
    }
    else
        panic("Error reading directory");
    return path_to_files;
}

// print details of input dataset
void print_dataset_details(vector<string> path_to_files, tTimer &timer)
{
    int average_contigs = 0;
    long int sequence_len;
    int num_contigs;

    for (int i = 0; i < path_to_files.size(); i++)
    {
        //Initial parse of file to get number of contigs and sequence length
        const char *filename = path_to_files[i].c_str();
        read_sequence(filename, sequence_len, num_contigs);
        printf("File Name: %s | Sequence Length: %ld | Number of Contigs: %d \n",
               (filename), sequence_len, num_contigs);
        average_contigs += num_contigs;
    }
    average_contigs = average_contigs / (int)path_to_files.size();

    cout << endl
         << "Average Contig Length: " << average_contigs << endl
         << endl;
    cout << "######################" << endl;
    timer.print("ALL FASTA READ");
    cout << "######################" << endl
         << endl;
}

// get contig block start and end indexes
vector<pair<int, int>> get_contig_blocks(int contig_size, int contig_threads, int kmer_len)
{
    int block_length = contig_size / contig_threads;
    int left_limit = 0, right_limit = block_length;
    vector<pair<int, int>> range_tuple;

    if (block_length < kmer_len)
    {
        range_tuple.push_back(make_pair(0, contig_size));
        return range_tuple;
    }

    range_tuple.push_back(make_pair(left_limit, right_limit));

    while ((right_limit - kmer_len + 1) + block_length <= contig_size and (right_limit - kmer_len + 1) > 0)
    {
        left_limit = right_limit - kmer_len + 1;
        right_limit = left_limit + block_length;
        range_tuple.push_back(make_pair(left_limit, right_limit));
    }
    if (contig_size - right_limit + kmer_len - 1 >= kmer_len)
    {
        left_limit = right_limit - kmer_len + 1;
        right_limit = contig_size;
        range_tuple.push_back(make_pair(left_limit, right_limit));
    }

    return range_tuple;
}

// Serial Implementation of DNA K-mer counting
void serial_kmer_counter(vector<string> path_to_files, tTimer &timer, int kmer_len)
{
    cout << "Serial Computation" << endl;
    for (int k = 0; k < path_to_files.size(); k++)
    {
        const char *filename = path_to_files[k].c_str();
        //Fasta reader part
        FASTAFILE *ffp = OpenFASTA(filename);
        char *sequence;
        char *name;
        int sequence_len;

        //create counter to hold kmers
        unordered_map<string, int> kmer_counter;

        while (ReadFASTA(ffp, &sequence, &name, &sequence_len))
        {
            // SAFE Cast!
            string t(sequence);

            for (int i = 0; i != sequence_len - kmer_len + 1; i++)
                kmer_counter[t.substr(i, kmer_len)]++;

            //for (auto &x : kmer_counter)
            //    std::cout << x.first << ": " << x.second << std::endl;

            free(sequence);
            free(name);
        }

        CloseFASTA(ffp);
    }
    timer.print("Completed in");
}

// openmp - nested - FASTA + Contig
// Level 1 - FASTA files are worked on in parallel
// Level 2 - Every contig within each FASTA file is distributed to 'contig_threads'
void openmp_kmer_nested(vector<string> path_to_files, tTimer &timer,
                        int kmer_len, int fasta_threads, int contig_threads)
{
    //Setting nested parallelism - Spawns (fasta_threads * contig_threads) threads
    omp_set_nested(1);
    cout << "OpenMP Nested Parallelism - FASTA (Level 1) + Contig (Level 2)" << endl;
#pragma omp parallel for num_threads(fasta_threads) default(shared) firstprivate(fasta_threads, contig_threads) schedule(dynamic, 1)
    for (int k = 0; k < path_to_files.size(); k++)
    {
        const char *filename = path_to_files[k].c_str();
        //Fasta reader part
        FASTAFILE *ffp = OpenFASTA(filename);
        char *sequence;
        char *name;
        int sequence_len;

        vector<string> contig_sequences;
        vector<int> contig_lengths;
        vector<string> contig_names;

        unordered_map<string, int> fasta_kmer_counter;

        while (ReadFASTA(ffp, &sequence, &name, &sequence_len))
        {
            //SAFE Cast!
            string s(sequence);
            contig_sequences.push_back(s);

            string n(name);
            contig_names.push_back(n);

            contig_lengths.push_back(sequence_len);

            free(sequence);
            free(name);
        }

        contig_threads = min(contig_threads, (int)contig_sequences.size());

#pragma omp parallel for num_threads(contig_threads) default(shared) schedule(dynamic, 1)
        for (int i = 0; i < contig_sequences.size(); i++)
        {
            for (int m = 0; m < contig_sequences[i].size() - kmer_len; m++)
                #pragma omp critical
                fasta_kmer_counter[contig_sequences[i].substr(m, kmer_len)]++;
        }
        CloseFASTA(ffp);
    }
    timer.print("Completed in");
}

// openmp - nested - FASTA + Contig Blocks
// Level 1 - FASTA files are worked on in parallel
// Level 2 - Every contig block within each FASTA file is distributed to 'contig_threads'
void openmp_kmer_nested_blocking(vector<string> path_to_files, tTimer &timer,
                                 int kmer_len, int fasta_threads, int contig_threads)
{
    //Setting nested parallelism - Spawns (fasta_threads * contig_threads) threads
    omp_set_nested(1);
    cout << "OpenMP Nested Parallelism - FASTA + Contig Blocking" << endl;
#pragma omp parallel for num_threads(fasta_threads) default(shared) firstprivate(contig_threads) schedule(dynamic, 1)
    for (int k = 0; k < path_to_files.size(); k++)
    {
        const char *filename = path_to_files[k].c_str();
        //Fasta reader part
        FASTAFILE *ffp = OpenFASTA(filename);
        char *sequence;
        char *name;
        int sequence_len;

        vector<string> contig_sequences;
        vector<int> contig_lengths;
        vector<string> contig_names;

        vector<unordered_map<string, int>> fasta_kmer_counter;

        while (ReadFASTA(ffp, &sequence, &name, &sequence_len))
        {
            //SAFE Cast!
            string s(sequence);
            contig_sequences.push_back(s);

            string n(name);
            contig_names.push_back(n);

            contig_lengths.push_back(sequence_len);

            free(sequence);
            free(name);
        }

        //contig_threads = min(contig_threads, (int)contig_sequences.size());
        //cout << filename << " | " \
             << "using " << contig_threads << " threads" << endl;

        for (int i = 0; i < contig_sequences.size(); i++)
        {
            //create counter to hold kmers
            unordered_map<string, int> contig_block_kmer_counter;

            vector<pair<int, int>> range_tuple =
                get_contig_blocks((int)contig_sequences[i].size(), contig_threads, kmer_len);

            vector<unordered_map<string, int>> contig_blocks_store(range_tuple.size());

            //#pragma omp declare reduction (merge : vector<unordered_map<string, int>> : \
                      omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

#pragma omp parallel for num_threads(contig_threads) default(shared) schedule(dynamic, 1) private(contig_block_kmer_counter)
            //reduction(merge: fasta_kmer_counter)
            for (int j = 0; j < range_tuple.size(); j++)
            {
                for (int m = range_tuple[j].first; m < range_tuple[j].second - kmer_len; m++)
                    contig_block_kmer_counter[contig_sequences[i].substr(m, kmer_len)]++;
                contig_blocks_store[j] = contig_block_kmer_counter;
            }
        }
        CloseFASTA(ffp);
    }
    timer.print("Completed in");
}

// openmp kmer computation - Every contig's kmers are computed in parallel
// but they are not aggregated in a final step
void openmp_kmer_computation_only(vector<string> path_to_files, tTimer &timer, int kmer_len, int contig_threads)
{
    cout << "OpenMP - Computation Only" << endl;
    for (int k = 0; k < path_to_files.size(); k++)
    {
        const char *filename = path_to_files[k].c_str();
        //Fasta reader part
        FASTAFILE *ffp = OpenFASTA(filename);
        char *sequence;
        char *name;
        int sequence_len;

        vector<string> contig_sequences;
        vector<int> contig_lengths;
        vector<string> contig_names;

        while (ReadFASTA(ffp, &sequence, &name, &sequence_len))
        {
            //SAFE Cast!
            string s(sequence);
            contig_sequences.push_back(s);

            string n(name);
            contig_names.push_back(n);

            contig_lengths.push_back(sequence_len);

            free(sequence);
            free(name);
        }

        //cout << filename << " | " \
             << "using " << contig_threads << " threads" << endl;

        for (int i = 0; i < contig_sequences.size(); i++)
        {
            unordered_map<string, int> contig_block_kmer_counter;

            vector<pair<int, int>> range_tuple =
                get_contig_blocks((int)contig_sequences[i].size(), contig_threads, kmer_len);

#pragma omp parallel for num_threads(contig_threads) default(shared) schedule(dynamic, 1) private(contig_block_kmer_counter)
            for (int j = 0; j < range_tuple.size(); j++)
                for (int m = range_tuple[j].first; m < range_tuple[j].second - kmer_len; m++)
                    contig_block_kmer_counter[contig_sequences[i].substr(m, kmer_len)]++;
        }
        CloseFASTA(ffp);
    }
    timer.print("Completed in");
}

// Driver function that calls kmer counter implementations with correct args
void driver(char *argv[])
{
    tTimer timer;
    int partitions = 4;
    int kmer_len = 17;
    int fasta_threads = 10;
    int contig_threads = 4;
    omp_set_dynamic(0);

    vector<string> path_to_files = read_directory(argv[1]);

    // print dataset contents
    print_dataset_details(path_to_files, timer);

    serial_kmer_counter(path_to_files, timer, kmer_len);
    openmp_kmer_nested(path_to_files, timer, kmer_len, fasta_threads, contig_threads);
    openmp_kmer_nested_blocking(path_to_files, timer, kmer_len, fasta_threads, contig_threads);
    openmp_kmer_computation_only(path_to_files, timer, kmer_len, contig_threads);
}

// single point of entry into kmer counting
int main(int argc, char *argv[])
{
    if (argc != 2)
        panic("missing data file name");

    driver(argv);

    return 0;
}
