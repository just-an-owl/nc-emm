

using System.Globalization;

namespace nc_emm
{
    internal class Score_Bins
    {
        public double[] bins_count;
        public double[] bins_boundaries;
        public Score_Bins(int bins)
        {
            this.bins_count = new double[bins + 1];
            this.bins_boundaries = new double[bins];
        }
        public void writeToFile(string filename)//tsv format
        {
            using (StreamWriter sw = new StreamWriter(filename))
            {
                int i = 0;
                sw.WriteLine("index\tCount\tBoundaries");
                for (i = 0; i < bins_count.Length-1; i++)
                {
                    sw.WriteLine(i.ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" + bins_count[i].ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" + bins_boundaries[i].ToString(CultureInfo.CreateSpecificCulture("en-GB")));
                }
              //  i += 1;
                sw.WriteLine(i.ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" + bins_count[i].ToString(CultureInfo.CreateSpecificCulture("en-GB")) + "\t" + bins_boundaries[i-1].ToString(CultureInfo.CreateSpecificCulture("en-GB")));//последняя граница дублируется в записи файла

            }
        }
    }
}
