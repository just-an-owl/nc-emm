

using System.Globalization;

namespace nc_emm
{
    internal class NucleotidChain
    {
        public string Name { get; }
        public string Chain { get; }
        private double[,] Frequency_matrix = new double[4, 4];
        private double[,] Matrix_triplet_frequencies = new double[16, 4];
        public double Euclidean_score { get; set; }
        public double Markov_score { get; set; }
        /// <summary>
        ///  A C T G
        /// A 
        /// C
        /// T
        /// G
        /// triplet
        ///   A C T G
        /// AA
        /// AC
        /// AT
        /// AG
        /// CA
        /// CC
        /// CT
        /// CG
        /// TA
        /// TC
        /// TT
        /// TG
        /// GA
        /// GC
        /// GT
        /// GG

        /// </summary>
        /// <param name="name"></param>
        /// <param name="chain"></param>
        /// <param name="frequency_matrix"></param>

        public NucleotidChain(string name, string chain, double[,] frequency_matrix)
        {
            Name = name;
            Chain = chain;
            Frequency_matrix = frequency_matrix;
        }
        public NucleotidChain(string name, string chain)
        {
            Name = name;
            Chain = chain;
            Markov_score = 0;
            Euclidean_score = 0;
        }
        public void create_frequency_matrix() {
            int[] counts = new int[4];//A C T G

            for (int i = 0; i < 4; i++)
            {
                counts[i] = 0;
                this.Frequency_matrix[0, i] = 0;
                this.Frequency_matrix[1, i] = 0;
                this.Frequency_matrix[2, i] = 0;
                this.Frequency_matrix[3, i] = 0;
            }
            for (int i = 0; i < Chain.Length-1; i++) {
                if (this.Chain[i] == 'A')
                {
                    counts[0] += 1;
                    if (this.Chain[i+1] == 'A')
                    {
                        this.Frequency_matrix[0, 0] += 1;
                    }
                    if (this.Chain[i + 1] == 'C')
                    {
                        this.Frequency_matrix[0, 1] += 1;
                    }
                    if (this.Chain[i + 1] == 'T')
                    {
                        this.Frequency_matrix[0, 2] += 1;
                    }
                    if (this.Chain[i + 1] == 'G')
                    {
                        this.Frequency_matrix[0, 3] += 1;
                    }
                }
                if (this.Chain[i] == 'C')
                {
                    counts[1] += 1;
                    if (this.Chain[i + 1] == 'A')
                    {
                        this.Frequency_matrix[1, 0] += 1;
                    }
                    if (this.Chain[i + 1] == 'C')
                    {
                        this.Frequency_matrix[1, 1] += 1;
                    }
                    if (this.Chain[i + 1] == 'T')
                    {
                        this.Frequency_matrix[1, 2] += 1;
                    }
                    if (this.Chain[i + 1] == 'G')
                    {
                        this.Frequency_matrix[1, 3] += 1;
                    }
                }
                if (this.Chain[i] == 'T')
                {
                    counts[2] += 1;
                    if (this.Chain[i + 1] == 'A')
                    {
                        this.Frequency_matrix[2, 0] += 1;
                    }
                    if (this.Chain[i + 1] == 'C')
                    {
                        this.Frequency_matrix[2, 1] += 1;
                    }
                    if (this.Chain[i + 1] == 'T')
                    {
                        this.Frequency_matrix[2, 2] += 1;
                    }
                    if (this.Chain[i + 1] == 'G')
                    {
                        this.Frequency_matrix[2, 3] += 1;
                    }
                }
                if (this.Chain[i] == 'G')
                {
                    counts[3] += 1;
                    if (this.Chain[i + 1] == 'A')
                    {
                        this.Frequency_matrix[3, 0] += 1;
                    }
                    if (this.Chain[i + 1] == 'C')
                    {
                        this.Frequency_matrix[3, 1] += 1;
                    }
                    if (this.Chain[i + 1] == 'T')
                    {
                        this.Frequency_matrix[3, 2] += 1;
                    }
                    if (this.Chain[i + 1] == 'G')
                    {
                        this.Frequency_matrix[3, 3] += 1;
                    }
                }
            }
            //досчитывание частоты для последнего нуклеотида в цепи
            if (this.Chain[this.Chain.Length - 1] == 'A') { 
                counts[0] += 1;
                this.Frequency_matrix[0, 0] += 1;
                this.Frequency_matrix[0, 1] += 1;
                this.Frequency_matrix[0, 2] += 1;
                this.Frequency_matrix[0, 3] += 1;
            }
            if (this.Chain[this.Chain.Length - 1] == 'C')
            {
                counts[1] += 1;
                this.Frequency_matrix[1, 0] += 1;
                this.Frequency_matrix[1, 1] += 1;
                this.Frequency_matrix[1, 2] += 1;
                this.Frequency_matrix[1, 3] += 1;
            }
            if (this.Chain[this.Chain.Length - 1] == 'T')
            {
                counts[2] += 1;
                this.Frequency_matrix[2, 0] += 1;
                this.Frequency_matrix[2, 1] += 1;
                this.Frequency_matrix[2, 2] += 1;
                this.Frequency_matrix[2, 3] += 1;
            }
            if (this.Chain[this.Chain.Length - 1] == 'G')
            {
                counts[3] += 1;
                this.Frequency_matrix[3, 0] += 1;
                this.Frequency_matrix[3, 1] += 1;
                this.Frequency_matrix[3, 2] += 1;
                this.Frequency_matrix[3, 3] += 1;
            }
            //расчет частоты
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++)
                {
                    if (counts[i] != 0)
                    {
                        this.Frequency_matrix[i, j] = this.Frequency_matrix[i, j] / counts[i];
                    }
                }
            }
        }
        public void create_triplet_frequency_matrix()
        {
            int[] counts = new int[16];

            for (int i = 0; i < 16; i++)
            {
                counts[i] = 0;
                this.Matrix_triplet_frequencies[i, 0] = 0;
                this.Matrix_triplet_frequencies[i, 1] = 0;
                this.Matrix_triplet_frequencies[i, 2] = 0;
                this.Matrix_triplet_frequencies[i, 3] = 0;
            }
            for (int i = 0; i < Chain.Length - 2; i++)
            {
                if (this.Chain[i].ToString() + this.Chain[i+1].ToString() == "AA")
                {
                    counts[0] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[0, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[0, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[0, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[0, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "AC")
                {
                    counts[1] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[1, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[1, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[1, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[1, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "AT")
                {
                    counts[2] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[2, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[2, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[2, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[2, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "AG")
                {
                    counts[3] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[3, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[3, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[3, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[3, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "CA")
                {
                    counts[4] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[4, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[4, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[4, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[4, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "CC")
                {
                    counts[5] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[5, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[5, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[5, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[5, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "CT")
                {
                    counts[6] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[6, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[6, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[6, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[6, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "CG")
                {
                    counts[7] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[7, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[7, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[7, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[7, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "TA")
                {
                    counts[8] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[8, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[8, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[8, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[8, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "TC")
                {
                    counts[9] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[9, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[9, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[9, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[9, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "TT")
                {
                    counts[10] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[10, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[10, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[10, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[10, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "TG")
                {
                    counts[11] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[11, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[11, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[11, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[11, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "GA")
                {
                    counts[12] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[12, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[12, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[12, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[12, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "GC")
                {
                    counts[13] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[13, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[13, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[13, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[13, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "GT")
                {
                    counts[14] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[14, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[14, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[14, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[14, 3] += 1;
                    }
                }
                if (this.Chain[i].ToString() + this.Chain[i + 1].ToString() == "GG")
                {
                    counts[15] += 1;
                    if (this.Chain[i + 2] == 'A')
                    {
                        this.Matrix_triplet_frequencies[15, 0] += 1;
                    }
                    if (this.Chain[i + 2] == 'C')
                    {
                        this.Matrix_triplet_frequencies[15, 1] += 1;
                    }
                    if (this.Chain[i + 2] == 'T')
                    {
                        this.Matrix_triplet_frequencies[15, 2] += 1;
                    }
                    if (this.Chain[i + 2] == 'G')
                    {
                        this.Matrix_triplet_frequencies[15, 3] += 1;
                    }
                }

            }


            //досчитывание частоты для последнего нуклеотида в цепи
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "AA")
            {
                counts[0] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[0, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "AC")
            {
                counts[1] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[1, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "AT")
            {
                counts[2] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[2, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "AG")
            {
                counts[3] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[3, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "CA")
            {
                counts[4] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[4, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "CC")
            {
                counts[5] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[5, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "CT")
            {
                counts[6] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[6, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "CG")
            {
                counts[7] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[7, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "TA")
            {
                counts[8] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[8, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "TC")
            {
                counts[9] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[9, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "TT")
            {
                counts[10] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[10, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "TG")
            {
                counts[11] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[11, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "GA")
            {
                counts[12] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[12, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "GC")
            {
                counts[13] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[13, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "GT")
            {
                counts[14] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[14, i] += 1;
            }
            if (this.Chain[this.Chain.Length - 2].ToString() + this.Chain[this.Chain.Length - 1].ToString() == "GG")
            {
                counts[15] += 1;
                for (int i = 0; i < 4; i++)
                    this.Matrix_triplet_frequencies[15, i] += 1;
            }

            //расчет частоты
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (counts[i] != 0)
                    {
                        this.Matrix_triplet_frequencies[i, j] = this.Matrix_triplet_frequencies[i, j] / counts[i];
                    }
                }
            }
        }
        public double[,] get_frequency_matrix()
        {
            return this.Frequency_matrix;
        }
        public string get_data_to_string()
        {
            string answer = "";

            answer += this.Name + "\t" 
                + this.Chain.Length.ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" 
                + Markov_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" 
                + Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB"));

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++) {
                    answer += "\t";
                    answer += this.Frequency_matrix[i, j].ToString(CultureInfo.CreateSpecificCulture("en-GB"));
                }
            }

            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    answer += "\t";
                    answer += this.Matrix_triplet_frequencies[i, j].ToString(CultureInfo.CreateSpecificCulture("en-GB"));
                }
            }

            return answer;
        }
    }
}
