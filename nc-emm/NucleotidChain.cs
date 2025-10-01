

namespace nc_emm
{
    internal class NucleotidChain
    {
        public string Name { get; }
        public string Chain { get; }
        private double[,] Frequency_matrix = new double[4, 4];
        public double Euclidean_score { get; set; }
        public double Markov_score { get; set; }
        /// <summary>
        ///  A C T G
        /// A 
        /// C
        /// T
        /// G
        /// 
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
        public double[,] get_frequency_matrix()
        {
            return this.Frequency_matrix;
        }
    }
}
