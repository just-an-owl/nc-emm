using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace nc_emm
{
    internal class ROC_score
    {
        public double sensitivity;
        public double specify;
        public double FPR;
        public int separation_threshold;
        public ROC_score(double specify, double sensetivity, int separation_threshold) {
            this.sensitivity = sensetivity;
            this.specify = specify; 
            this.separation_threshold = separation_threshold;
            this.FPR = 1 - specify;
        }

    }
}
