using System;
using System.Collections.Generic;
using System.Globalization;
using System.Diagnostics;
using System.IO;





namespace nc_emm
{
    internal class Matrix_interactions
    {
        private List<NucleotidChain> _referens;
        private List<NucleotidChain> _testable_circuits;
        private List<List<double>> _matrix_difference;
        private double[,] _mean_matrix;
        private double[,] _standard_deviation;
        private Score_Bins referens_Bins;
        private Score_Bins testable_Bins;
        private int point_of_separation;
        private double F1_score;

        private double max_enh = 0;
        
        private double max_prom = 0;
        public Matrix_interactions(List<NucleotidChain> nucleotidChains) => _referens = new List<NucleotidChain>(nucleotidChains);
        public Matrix_interactions() { }
        public Matrix_interactions(string referens_path, string testeble_path, int bins)//готовый пайплайн работы проекта
        {
            Stopwatch stopwatch = new Stopwatch();

            Console.WriteLine("Analyzing the reference file");
            stopwatch.Start();
            
            this.read_from_file(referens_path, 0);
            
            stopwatch.Stop();
            TimeSpan timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            stopwatch = new Stopwatch();
            Console.WriteLine("Analyzing the testeble file");
            stopwatch.Start();

            this.read_from_file(testeble_path, 1);

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            
            Console.WriteLine("Creating a model matrix");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            this.create_mean_matrix();
            
            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            stopwatch = new Stopwatch();

            Console.WriteLine("Creating a matrix of standard deviations");
            stopwatch = new Stopwatch();
            stopwatch.Start();
            
            this.calculate_standard_deviation();
            
            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            Console.WriteLine("Matrix of standart deviation:");
            string Matrix = "";
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Matrix += this._standard_deviation[i, j];
                    Matrix += " ";
                }
                Matrix += "\n";
            }
            Console.WriteLine(Matrix);

            Console.WriteLine("Calculating the Euclidean distance");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            this.Euclidean_distance_calculation(0);

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            this.Euclidean_distance_calculation(1);

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");

            Console.WriteLine("Find change dot");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            this.find_change_dot(bins);

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine("Chenge dot is " + this.point_of_separation.ToString());
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");


            Console.WriteLine("Calculating F1 score");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            this.F1();

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");

            Console.WriteLine("Calculating ROC");
            stopwatch = new Stopwatch();
            stopwatch.Start();

            List<ROC_score> rOC_Scores = this.ROC();

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");
            
            Console.WriteLine("Writing files");
            stopwatch = new Stopwatch();
            stopwatch.Start();
            
            DateTime dateTime = DateTime.Now;
            string dir_path = dateTime.ToString("dd.MM.yyyy hh.mm.ss");
            Directory.CreateDirectory(dir_path);

            this.writeToFilePromoteScore(dir_path+ "/TestableCircuitsScore.txt");
            this.writeToFileEnhancerScore(dir_path + "/ReferensScore.txt");
            this.writeToFileStatisticalData(dir_path + "/StatisticalData.txt");
            this.writeToFileROCData(dir_path + "/ROCData.tsv", rOC_Scores);
            this.writeMatrix(dir_path+ "/ModelData.tsv", this._mean_matrix);
            this.writeMatrix(dir_path + "/StandartDeviationData.tsv", this._standard_deviation);
            this.referens_Bins.writeToFile(dir_path + "/ReferensBins.tsv");
            this.testable_Bins.writeToFile(dir_path + "/TestableBins.tsv");

            stopwatch.Stop();
            timeSpan = stopwatch.Elapsed;
            Console.WriteLine($"Lead time: {timeSpan.TotalSeconds} секунд");

        }
        public void create_matrix_differences()
        {
            this._matrix_difference = new List<List<double>>();
            for (int i = 0; i < _referens.Count-1; i++) {
                _matrix_difference.Add(new List<double>());
                for (int j = i + 1; j < _referens.Count; j++) {
                    _matrix_difference[i].Add(counting_mean_matrix_difference(_referens[i].get_frequency_matrix(), _referens[j].get_frequency_matrix()));
                }
                if (i%10_000 == 0)
                {
                    System.Console.WriteLine($"Расчитано {i} / {_referens.Count} первичных матриц");
                }
            }
        }
        private double counting_mean_matrix_difference(double[,] f_matrix,  double[,] s_matrix)
        {
            double[] del_matrix = new double[4*4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    del_matrix[4*i+j] = f_matrix[i, j] - s_matrix[i, j];
                }
            }
            double s = del_matrix.Sum();
            return s / del_matrix.Length;
        }
        public void read_from_file(string filename, int mod)
        {
            if (mod == 0)
            {
                this._referens = new List<NucleotidChain>();
                using (StreamReader reader = new StreamReader(filename))
                {
                    while (reader.EndOfStream != true)
                    {
                        this._referens.Add(new NucleotidChain(reader.ReadLine(), reader.ReadLine().ToUpper()));
                    }
                }
                this._referens.ForEach(x => { x.create_frequency_matrix(); });
            }
            if (mod == 1)
            {
                this._testable_circuits = new List<NucleotidChain>();
                using (StreamReader reader = new StreamReader(filename))
                {
                    while (reader.EndOfStream != true)
                    {
                        this._testable_circuits.Add(new NucleotidChain(reader.ReadLine(), reader.ReadLine().ToUpper()));
                    }
                }
                this._testable_circuits.ForEach(x => { x.create_frequency_matrix(); });
            }
        }
        public void create_mean_matrix()
        {
            this._mean_matrix = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                this._mean_matrix[0, i] = 0;
                this._mean_matrix[1, i] = 0;
                this._mean_matrix[2, i] = 0;
                this._mean_matrix[3, i] = 0;
            }
            for (int i = 0; i < this._referens.Count; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        this._mean_matrix[j,k] += this._referens[i].get_frequency_matrix()[j, k];
                    }
                }
            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    this._mean_matrix[j, k] /= this._referens.Count;
                }
            }

        }
        public double[,] get_mean_matrix() { return this._mean_matrix; }
        public void calculate_standard_deviation()//среднеквадратичное отклонение
        {
            this._standard_deviation = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                this._standard_deviation[0, i] = 0;
                this._standard_deviation[1, i] = 0;
                this._standard_deviation[2, i] = 0;
                this._standard_deviation[3, i] = 0;
            }
            for (int i = 0; i < this._referens.Count; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        this._standard_deviation[j, k] += System.Math.Pow(this._mean_matrix[j,k] - this._referens[i].get_frequency_matrix()[j, k],2);
                    }
                }

            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    this._standard_deviation[j, k] /= this._referens.Count;
                }
            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    this._standard_deviation[j, k] = System.Math.Pow(this._standard_deviation[j,k], 0.5);
                }
            }
        }
        public double[,] get_standard_deviation() {  return this._standard_deviation; }
        public double[,] get_norm_standard_deviation()
        {
            double[,] norm_standard_ = new double[4, 4];

            for (int i = 0; i < 4; i++)
            {
                norm_standard_[0, i] = 0;
                norm_standard_[1, i] = 0;
                norm_standard_[2, i] = 0;
                norm_standard_[3, i] = 0;
            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    norm_standard_[j, k] = this._standard_deviation[j,k] /this._mean_matrix[j, k];
                }
            }
            return norm_standard_;
        }

        /*
        public void writeToFileEnhancerScore(string path, int min_length, float min_score)
        {
            // полная перезапись файла 
            using (StreamWriter writer = new StreamWriter(path, false))
            {
                for (int i = 0; i < this._referens.Count; i++)
                    if (this._testable_circuits[i].Euclidean_score >= min_score && this._referens[i].Chain.Length >= min_length)
                        writer.WriteLine(this._referens[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }*/
        
        /*public async void writeToFilePromoteScore(string path, Boolean isExthist)
        {
            using (StreamWriter writer = new StreamWriter(path, true))
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                    await writer.WriteLineAsync(this._testable_circuits[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }*/
        public void Euclidean_distance_calculation(int mod)
        {
            if (mod == 0)
            {
                for (int i = 0; i < this._referens.Count; i++)
                {
                    double calc_score = 0;
                    for (int j = 0; j < 4; j++)
                    {
                        for (int k = 0; k < 4; k++)
                        {

                            calc_score += Math.Pow((this._referens[i].get_frequency_matrix()[j, k]) - (this._mean_matrix[j, k]), 2);

                        }

                    }
                    this._referens[i].Euclidean_score = Math.Sqrt(calc_score);
                    if (this.max_enh < this._referens[i].Euclidean_score)
                        this.max_enh = this._referens[i].Euclidean_score;
                }
            }
            if (mod == 1)
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    double calc_score = 0;
                    for (int j = 0; j < 4; j++)
                    {
                        for (int k = 0; k < 4; k++)
                        {

                            calc_score += Math.Pow((this._testable_circuits[i].get_frequency_matrix()[j, k]) - (this._mean_matrix[j, k]), 2);

                        }

                    }
                    this._testable_circuits[i].Euclidean_score = Math.Sqrt(calc_score);
                    if (this.max_prom < this._testable_circuits[i].Euclidean_score)
                        this.max_prom = this._testable_circuits[i].Euclidean_score;
                }
            }
            }
        public void Log_Euclidean_distance_calculation(int mod)
        {
            if (mod == 0) { 
                for (int i = 0; i < this._referens.Count; i++)
                {
                    double calc_score = 0;
                    for (int j = 0; j < 4; j++)
                    {
                        for (int k = 0; k < 4; k++) {
                            if (this._referens[i].get_frequency_matrix()[j, k] != 0)
                                calc_score += Math.Pow(Math.Log(this._referens[i].get_frequency_matrix()[j, k]) - Math.Log(this._mean_matrix[j, k]), 2);
                            else
                                calc_score += Math.Pow(Math.Log(this._mean_matrix[j, k]), 2);
                        }

                    }
                    this._referens[i].Euclidean_score = Math.Sqrt(calc_score);
                }
            }
        }
        public void Counting_Markov_distance() {
            for (int i = 0; i < this._referens.Count; i++) {
                double score = 1;
                for (int j = 1; j < this._referens[i].Chain.Length; j++) {
                    if (this._referens[i].Chain[j-1] == 'A')
                    {
                        if (this._referens[i].Chain[j] == 'A')
                        {
                            score *= Math.Log(this._mean_matrix[0, 0]);
                        }
                        if (this._referens[i].Chain[j] == 'C')
                        {
                            score *= Math.Log(this._mean_matrix[0, 1]);
                        }
                        if (this._referens[i].Chain[j] == 'T')
                        {
                            score *= Math.Log(this._mean_matrix[0, 2]);
                        }
                        if (this._referens[i].Chain[j] == 'G')
                        {
                            score *= Math.Log(this._mean_matrix[0, 3]);
                        }
                    }
                    if (this._referens[i].Chain[j - 1] == 'C')
                    {
                        if (this._referens[i].Chain[j] == 'A')
                        {
                            score *= Math.Log(this._mean_matrix[1, 0]);
                        }
                        if (this._referens[i].Chain[j] == 'C')
                        {
                            score *= Math.Log(this._mean_matrix[1, 1]);
                        }
                        if (this._referens[i].Chain[j] == 'T')
                        {
                            score *= Math.Log(this._mean_matrix[1, 2]);
                        }
                        if (this._referens[i].Chain[j] == 'G')
                        {
                            score *= Math.Log(this._mean_matrix[1, 3]);
                        }
                    }
                    if (this._referens[i].Chain[j - 1] == 'T')
                    {
                        if (this._referens[i].Chain[j] == 'A')
                        {
                            score *= Math.Log(this._mean_matrix[2, 0]);
                        }
                        if (this._referens[i].Chain[j] == 'C')
                        {
                            score *= Math.Log(this._mean_matrix[2, 1]);
                        }
                        if (this._referens[i].Chain[j] == 'T')
                        {
                            score *= Math.Log(this._mean_matrix[2, 2]);

                        }
                        if (this._referens[i].Chain[j] == 'G')
                        {
                            score *= Math.Log(this._mean_matrix[2, 3]);

                        }
                    }
                    if (this._referens[i].Chain[j - 1] == 'G')
                    {
                        if (this._referens[i].Chain[j] == 'A')
                        {
                            score *= Math.Log(this._mean_matrix[3, 0]);

                        }
                        if (this._referens[i].Chain[j] == 'C')
                        {
                            score *= Math.Log(this._mean_matrix[3, 1]);

                        }
                        if (this._referens[i].Chain[j] == 'T')
                        {
                            score *= Math.Log(this._mean_matrix[3, 2]);

                        }
                        if (this._referens[i].Chain[j] == 'G')
                        {
                            score *= Math.Log(this._mean_matrix[3, 3]);

                        }
                    }
                }
                this._referens[i].Markov_score = score * 1000 / this._referens[i].Chain.Length;
            }
        }
        public double get_mean_lenth_of_chain()
        {
            double score = 0;
            for (int i = 0; i < this._referens.Count; i++)
            {
                score += this._referens[i].Chain.Length;
            }
            return score / this._referens.Count;
        }
        public void sequence_consistency_calculation(string path)
        {
            List<double> in_score = new List<double>();
            List<double> out_score = new List<double>();
            for (int i = 0; i<this._referens.Count; i++)
            {

            }
        }
        public void get_score_sequences_with_length(string path, int minL, int maxL, int mode) {
            if (mode == 0)
            using (StreamWriter writer = new StreamWriter(path, false))
            {
                for (int i = 0; i < this._referens.Count; i++)
                    if (this._referens[i].Chain.Length>minL && this._referens[i].Chain.Length < maxL)
                        writer.WriteLine(this._referens[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }
        public int getCountPromoter() { return this._testable_circuits.Count;}
        
        private int checkBound(double[] boundaries, double check_num)
        {
            int position = 0;
            for (int i = 0; i < boundaries.Length; i++)
            {
                if (check_num < boundaries[i])
                {
                    return position;
                }
                position++;
            }
            return position;
        }
        private Score_Bins createQuntilBins(int mod, int bins)
        {
            //создание общих границ для дальнейшего разбиения
            double[] bins_count_procent = new double[bins + 1];
            double[] bins_boundaries = new double[bins];
            double step = Math.Max(this.max_enh, this.max_prom) / bins;
            bins_count_procent[0] = 0;
            for (int i = 0; i < bins_boundaries.Length; i++)
            {
                bins_boundaries[i] = step * (i + 1);
                bins_count_procent[i + 1] = 0;
            }
            if (mod == 0)
            {
                for (int i = 0; i < this._referens.Count; i++)
                {
                    bins_count_procent[checkBound(bins_boundaries, this._referens[i].Euclidean_score)] += 1;
                }
                for (int i = 0; i < bins_count_procent.Length; i++)//результаты преобразования приводятся к долевому виду
                {
                    bins_count_procent[i] = bins_count_procent[i] / this._referens.Count;
                }
                Score_Bins referens_bins = new Score_Bins(bins);
                referens_bins.bins_count = bins_count_procent;
                referens_bins.bins_boundaries = bins_boundaries;
                return referens_bins;
            }
            if (mod == 1)
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    bins_count_procent[checkBound(bins_boundaries, this._testable_circuits[i].Euclidean_score)] += 1;
                }
                for (int i = 0; i < bins_count_procent.Length; i++)//результаты преобразования приводятся к долевому виду
                {
                    bins_count_procent[i] = bins_count_procent[i] / this._testable_circuits.Count;
                }
                Score_Bins test_bins = new Score_Bins(bins);
                test_bins.bins_count = bins_count_procent;
                test_bins.bins_boundaries = bins_boundaries;
                return test_bins;
            }
            return null;
        }
        public double F1() {
            double precision = 0;
            double recall = 0;
            double tp = 0;
            double fp = 0;
            int refer_max_index = Array.IndexOf(this.referens_Bins.bins_count, this.referens_Bins.bins_count.Max());
            int testeble_max_index = Array.IndexOf(this.testable_Bins.bins_count, this.testable_Bins.bins_count.Max());
            if (refer_max_index > testeble_max_index)
            {
                for (int i = 0; i < this.point_of_separation; i++)
                {
                    tp += this.referens_Bins.bins_count[i];
                    fp += this.testable_Bins.bins_count[i];
                }
            }
                
            precision = tp / (tp + fp);
            recall = tp;

            this.F1_score = 2 * (precision * recall) / (precision + recall);

            return this.F1_score; }
        private List<double> calc_ROC_score(int point)
        {

            double TN = 0;
            double FP = 0;
            double TP = 0;
            double FN = 0;

            List<double> Roc_score = new List<double>();//Se,Sp

            for (int i = 0; i < this.referens_Bins.bins_count.Length; i++)
            {
                if (i <= point)
                {
                    TP += this.referens_Bins.bins_count[i];
                    FP += this.testable_Bins.bins_count[i];

                }
                if (i > point)
                {
                    FN += this.referens_Bins.bins_count[i];
                    TN += this.testable_Bins.bins_count[i];
                }

            }
            
            
            Roc_score.Add(TP / (TP + FN));
            Roc_score.Add(TN / (TN + FP));

            return Roc_score;
        }
        public List<ROC_score> ROC() {


            List<ROC_score> Roc_scores = new List<ROC_score>();
            for (int i = 0; i < this.referens_Bins.bins_count.Length; i++)
            {
                int separation_threshold = i;
                List<double> roc_score = calc_ROC_score(separation_threshold);
                Roc_scores.Add(new ROC_score(roc_score[1], roc_score[0],separation_threshold));
            }
            Roc_scores.Sort((p1, p2) => p1.FPR.CompareTo(p2.FPR));
            
            return Roc_scores; }//передает отсортированный по специфичности массив лист roc объектов
        public double find_change_dot(int bins) {

            Score_Bins refer_bins = this.createQuntilBins(0, bins);
            Score_Bins tetstable_Bins = this.createQuntilBins(1, bins);
            this.referens_Bins = refer_bins;
            this.testable_Bins = tetstable_Bins;
            int dot = -1;
            int refer_max_index = Array.IndexOf(refer_bins.bins_count,refer_bins.bins_count.Max());
            int testeble_max_index = Array.IndexOf(tetstable_Bins.bins_count, tetstable_Bins.bins_count.Max());
            if(refer_max_index > testeble_max_index)
            {
                int i = 0;
                while (dot == -1)
                {
                    if (refer_bins.bins_count[i] - tetstable_Bins.bins_count[i] < 0)
                    {
                        dot = i;
                    }
                    
                    i++;

                }
            }
            else
            {
                int i = 0;
                while (dot == -1 && i<refer_bins.bins_count.Length)
                {
                    if (refer_bins.bins_count[i] - tetstable_Bins.bins_count[i] < 0)
                    {
                        dot = i;
                    }

                    i++;

                }
            }

            this.point_of_separation = dot;
            return dot; 
        }//поиск точки разделимости последовательностей



        //функции записи
        public void writeToFileEnhancerScore(string path)
        {
            // полная перезапись файла 
            using (StreamWriter writer = new StreamWriter(path, false))
            {
                for (int i = 0; i < this._referens.Count; i++)
                    writer.WriteLine(this._referens[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }
        public void writeToFilePromoteScore(string path)
        {
            // полная перезапись файла 
            using (StreamWriter writer = new StreamWriter(path, false))
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                    writer.WriteLine(this._testable_circuits[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }//запись всех полученых оценок
        public void writeToFilePromoteScore(string path, int min_length, float min_score)
        {
            // полная перезапись файла 
            using (StreamWriter writer = new StreamWriter(path, false))
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                    if (this._testable_circuits[i].Euclidean_score >= min_score && this._testable_circuits[i].Chain.Length >= min_length)
                        writer.WriteLine(this._testable_circuits[i].Euclidean_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        }//запись оценок от min_score с цепями от длины min_lenth
        public void writeLenthReferens(string path)
        {
            using (StreamWriter sw = new StreamWriter(path))
            {
                for (int i = 0; i < this._referens.Count; i++)
                {
                    sw.WriteLine(this._referens[i].Chain.Length);
                }
            }
        }
        public void writeLenthTesteble(string path)
        {
            using (StreamWriter sw = new StreamWriter(path))
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    sw.WriteLine(this._testable_circuits[i].Chain.Length);
                }
            }
        }
        public void writeQuntilfile(string path, int bins, int mod)
        {
            if (mod == 0)
            {
                int[] bins_count = new int[bins + 1];
                double[] bins_boundaries = new double[bins];
                double step = this.max_enh / bins;
                bins_count[0] = 0;
                for (int i = 0; i < bins_boundaries.Length; i++)
                {
                    bins_boundaries[i] = step * (i + 1);
                    bins_count[i + 1] = 0;
                }
                for (int i = 0; i < this._referens.Count; i++)
                {
                    bins_count[checkBound(bins_boundaries, this._referens[i].Euclidean_score)] += 1;
                }
                using (StreamWriter sw = new StreamWriter(path))
                {
                    for (int i = 0; i < bins_count.Length; i++)
                    {
                        //Console.WriteLine($"{((double)bins_count[i] / this._referens.Count * 100).ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                        sw.WriteLine($"{((double)bins_count[i] / this._referens.Count * 100).ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                    }
                }
            }
            if (mod == 1)
            {
                int[] bins_count = new int[bins + 1];
                double[] bins_boundaries = new double[bins];
                double step = this.max_prom / bins;
                bins_count[0] = 0;
                for (int i = 0; i < bins_boundaries.Length; i++)
                {
                    bins_boundaries[i] = step * (i + 1);
                    bins_count[i + 1] = 0;
                }
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    bins_count[checkBound(bins_boundaries, this._testable_circuits[i].Euclidean_score)] += 1;
                }
                using (StreamWriter sw = new StreamWriter(path))
                {
                    for (int i = 0; i < bins_count.Length; i++)
                        sw.WriteLine($"{((double)bins_count[i] / this._testable_circuits.Count * 100).ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                }
            }
            if (mod == 2)
            {
                int[] bins_count = new int[bins + 1];
                double[] bins_boundaries = new double[bins];
                double step = this.max_enh / bins;
                bins_count[0] = 0;

                for (int i = 0; i < bins_boundaries.Length; i++)
                {
                    bins_boundaries[i] = step * (i + 1);
                    bins_count[i + 1] = 0;
                }
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    bins_count[checkBound(bins_boundaries, this._testable_circuits[i].Euclidean_score)] += 1;
                }
                using (StreamWriter sw = new StreamWriter(path))
                {
                    for (int i = 0; i < bins_count.Length; i++)
                        sw.WriteLine($"{((double)bins_count[i] / this._testable_circuits.Count * 100).ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                }
                using (StreamWriter sw = new StreamWriter("boundaries" + path))
                {
                    for (int i = 0; i < bins_boundaries.Length; i++)
                        sw.WriteLine($"{bins_boundaries[i].ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                    sw.WriteLine($"{max_prom.ToString(CultureInfo.CreateSpecificCulture("en-GB"))}");
                }
            }
        }
        public void writeGCfrequency(string path, int mod)
        {
            using (StreamWriter sw = new StreamWriter(path))
            {
                if (mod == 0)
                {
                    for (int i = 0; i < this._referens.Count; i++)
                    {
                        sw.WriteLine(this._referens[i].get_frequency_matrix()[1, 3].ToString(CultureInfo.CreateSpecificCulture("en-GB")));
                    }
                }
                if (mod == 1)
                {
                    for (int i = 0; i < this._referens.Count; i++)
                    {
                        if (this._referens[i].Chain.Length > 150)
                            sw.WriteLine(this._referens[i].get_frequency_matrix()[1, 3].ToString(CultureInfo.CreateSpecificCulture("en-GB")));
                    }
                }
                if (mod == 2)
                {
                    for (int i = 0; i < this._testable_circuits.Count; i++)
                    {
                        if (this._testable_circuits[i].Chain.Length > 150)
                            sw.WriteLine(this._testable_circuits[i].get_frequency_matrix()[1, 3].ToString(CultureInfo.CreateSpecificCulture("en-GB")));
                    }
                }
            }
        }
        public void writeHightScorePromoter(string path, int min_length, double referens_score)
        {
            using (StreamWriter writer = new StreamWriter(path))
            {
                for (int i = 0; i < this._testable_circuits.Count; i++)
                {
                    if (this._testable_circuits[i].Euclidean_score >= referens_score && this._testable_circuits[i].Chain.Length >= min_length)
                    {
                        writer.WriteLine(this._testable_circuits[i].Name);
                        writer.WriteLine(this._testable_circuits[i].Chain);
                    }
                }
            }
        }
        public void createFileForCorrelationCheck(string path)
        {
            using (StreamWriter writer = new StreamWriter(path))
            {

            }
        }//не написана
        public void writeToFileStatisticalData(string path) {

            using (StreamWriter writer = new StreamWriter(path)) {
                writer.WriteLine("F1 score:");
                writer.WriteLine(this.F1_score.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
                writer.WriteLine(point_of_separation.ToString(CultureInfo.CreateSpecificCulture("en-GB")));
            }
        
        }
        public void writeToFileROCData(string path, List<ROC_score> ROC_scores)
        {
            using (StreamWriter writer= new StreamWriter(path))
            {
                writer.WriteLine("Index\tFPR\tSensitivity\tSpecify");
                for (int i = 0; i < ROC_scores.Count; i++) {
                    writer.WriteLine(i.ToString() + "\t" + ROC_scores[i].FPR.ToString(CultureInfo.CreateSpecificCulture("en-GB"))
                        + "\t" + ROC_scores[i].sensitivity.ToString(CultureInfo.CreateSpecificCulture("en-GB")) 
                        + "\t" + ROC_scores[i].specify.ToString(CultureInfo.CreateSpecificCulture("en-GB")));

                }
            }
        }
        public void writeMatrix(string path, double[,] matrix)//для квадратных матриц 4x4
        {
            string Matrix = "";
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Matrix += matrix[i, j];
                    Matrix += "\t";
                }
                Matrix += "\n";
            }
            using(StreamWriter writer= new StreamWriter(path)) {  writer.WriteLine(Matrix); }
        }
    }
}
