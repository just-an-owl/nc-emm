using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace nc_emm
{
    internal class CombineNonCodChain
    {
        private double[,] _mean_matrix;
        private List<NucleotidChain> _referens;

        public CombineNonCodChain() { }
        public void readAndWrite(string inputFile, string outputFile, int blockSize = 6000)
        {
            using(StreamWriter writer = new StreamWriter(outputFile,append:true))
            {
                writer.WriteLine("Name\tLength\tMarkov_score\tEuclidean_score\tAA\tAC\tAT\tAG\tCA\tCC\tCT\tCG\tTA\tTC\tTT\tTG\tGA\tGC\tGT\tGG" +
                            "\tAAA\tAAC\tAAT\tAAG\tACA\tACC\tACT\tACG\tATA\tATC\tATT\tATG\tAGA\tAGC\tAGT\tAGG\tCAA\tCAC\tCAT\tCAG\tCCA\tCCC\tCCT\tCCG\tCTA\tCTC\tCTT\tCTG\tCGA\tCGC\tCGT\tCGG\tTAA\tTAC\tTAT\tTAG\tTCA\tTCC\tTCT\tTCG\tTTA\tTTC\tTTT\tTTG\tTGA\tTGC\tTGT\tTGG\tGAA\tGAC\tGAT\tGAG\tGCA\tGCC\tGCT\tGCG\tGTA\tGTC\tGTT\tGTG\tGGA\tGGC\tGGT\tGGG");
            }
            //readModel("ModelData.tsv");
            string global_name = "";
            using (StreamReader reader = new StreamReader(inputFile))
            {
                while (reader.EndOfStream != true)
                {
                    global_name = reader.ReadLine();
                    string source = reader.ReadLine().ToUpper();
                    int startPos = 0;
                    while (startPos+blockSize <= source.Length)
                    {
                        string sub = source.Substring(startPos, blockSize);

                        NucleotidChain nucleotidChain = new NucleotidChain(global_name + "+" + startPos.ToString(), sub);
                        nucleotidChain.Euclidean_score = 0;
                        nucleotidChain.Markov_score = 0;
                        nucleotidChain.create_frequency_matrix();
                        nucleotidChain.create_triplet_frequency_matrix();
                        using (StreamWriter writer = new StreamWriter(outputFile,append:true))
                        {
                            writer.WriteLine(nucleotidChain.get_data_to_string(true));
                        }

                        startPos = startPos + 1000;
                    }
                
                }
            }

        }
        
        public void readModel(string inputFile)
        {
            using (StreamReader reader = new StreamReader(inputFile))
            {
                int n = 0;
                while (reader.EndOfStream != true)
                {
                    string[] cloude = reader.ReadLine().Split('\t');
                    for (int i = 0; i < cloude.Length; i++)
                    {
                        _mean_matrix[n,i] = double.Parse(cloude[i]);
                    }
                    n++;
                }
            }
        }
    }
}
