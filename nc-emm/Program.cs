using nc_emm;
using System;


namespace nc_emm
{
    class Program
    {
        static void Main(string[] args)
        {
            string testeble = null;
            string referens = null;
            bool triplet_mode = false;
            int bins = -1;
            Matrix_interactions matrix;
            // Проверка наличия аргументов
            if (args.Length == 0)
            {
                CombineNonCodChain combine = new CombineNonCodChain();
                combine.readAndWrite("seq_1kk_6000_random_2_-erna.txt", "background_1kk_6000_random.tsv");
                /*
                Console.WriteLine("Используется режжим тестового запуска на стандартных параметрах. Используйте -h чтобы увидеть справку о параметрах работы программы.");
                matrix = new Matrix_interactions("seq.txt", "seq_500k_700_random_2.txt", 10000, true);
                return;
                */
            }

            for (int i = 0; i < args.Length; i++)
            {
                switch (args[i])
                {
                    case "-h":
                    case "--help":
                        ShowHelp();
                        break;
                    case "-v":
                    case "--version":
                        ShowVersion();
                        break;
                    case "-r":
                    case "--referens":
                        if (i + 1 < args.Length)
                        {
                            referens = args[++i];
                        }
                        break;
                    case "-t":
                    case "--testeble":
                        if (i + 1 < args.Length)
                        {
                            testeble = args[++i];
                        }
                        break;
                    case "-b":
                    case "--bins":
                        if (i + 1 < args.Length)
                        {
                            bins = Convert.ToInt32(args[++i]);
                        }
                        break;
                    case "--triplet":
                        if (i + 1 < args.Length)
                        {
                            triplet_mode = Convert.ToBoolean(args[++i]);
                        }
                        break;
                    default:
                        Console.WriteLine($"Неизвестный параметр: {args[i]}");
                        break;
                }
            }
            if (referens != null && testeble != null && bins > 0)
            {
                matrix = new Matrix_interactions(referens, testeble, bins, triplet_mode);
            }
            else
            {
                Console.WriteLine("Пожалуйста, проверьте корректность переданных параметров");
            }
        }

        static void ShowHelp()
        {
            Console.WriteLine("Справка по программе:");
            Console.WriteLine("  -h, --help     Показать эту справку");
            Console.WriteLine("  -v, --version  Показать версию программы");
            Console.WriteLine("  -r, --referens Путь к файлу референсу");
            Console.WriteLine("  -t, --testeble Путь к файлу для тестирования");
            Console.WriteLine("  -b, --bins Количество разделений при кластеризации значений Евклидова расстояния");

        }

        static void ShowVersion()
        {
            Console.WriteLine("Версия 1.0.0");
        }
    }
}
