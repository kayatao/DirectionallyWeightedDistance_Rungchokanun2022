using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//For Stopwatch
using System.Diagnostics;
//For PointF
using System.Drawing;
//For IO
using System.IO;
//For Image Utilities
using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.CvEnum;
using Emgu.CV.UI;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class LAS_Liu2006
    {
        public Image<Gray, Single> Gx;
        public Image<Gray, Single> Gy;
        public Image<Gray, Single> Gxx;
        public Image<Gray, Single> Gyy;
        public Image<Gray, Single> Gxy;
        Image<Gray, Single> Theta;
        Image<Gray, Single> Coh;

        public Image<Gray, Single> POFx;
        public Image<Gray, Single> POFy;

        Dictionary<Point, Dictionary<string, Point[]>> observedPDirectionSector_ipj;
        Dictionary<Point, Dictionary<string, Point[]>> observedPDirectionSector_imj_1;


        public double sigmaAvgWindows = 5;
        public double sigmaGaussian = 20;
        public double radius = 60;
        public int numOfSectors = 50;
        public int blockSize = 4;

        public LAS_Liu2006()
        {

        }

        ~LAS_Liu2006()
        {

        }

        private void Initialize(Image<Gray, Single> input)
        {
            Gx = new Image<Gray, Single>(input.Width, input.Height);
            Gy = new Image<Gray, Single>(input.Width, input.Height);
            Gxx = new Image<Gray, Single>(input.Width, input.Height);
            Gyy = new Image<Gray, Single>(input.Width, input.Height);
            Gxy = new Image<Gray, Single>(input.Width, input.Height);
            Theta = new Image<Gray, Single>(input.Width, input.Height);
            Coh = new Image<Gray, Single>(input.Width, input.Height);


            POFx = new Image<Gray, Single>(input.Width, input.Height);
            POFy = new Image<Gray, Single>(input.Width, input.Height);

            if (observedPDirectionSector_ipj == null)
            {
                observedPDirectionSector_ipj = new Dictionary<Point, Dictionary<string, Point[]>>();
                observedPDirectionSector_imj_1 = new Dictionary<Point, Dictionary<string, Point[]>>();
            }

        }

        public void free()
        {
            Gx.Dispose();
            Gy.Dispose();
            Gxx.Dispose();
            Gyy.Dispose();
            Gxy.Dispose();
            Theta.Dispose();
            Coh.Dispose();

            POFx.Dispose();
            POFy.Dispose();

            foreach (KeyValuePair<Point, Dictionary<string, Point[]>> kvp in observedPDirectionSector_ipj)
            {
                observedPDirectionSector_ipj[kvp.Key].Clear();
            }
            foreach (KeyValuePair<Point, Dictionary<string, Point[]>> kvp in observedPDirectionSector_imj_1)
            {
                observedPDirectionSector_imj_1[kvp.Key].Clear();
            }

            observedPDirectionSector_ipj.Clear();
            observedPDirectionSector_imj_1.Clear();
        }

        public void runLAS(string fpImgPath, string savePath)
        {
            Stopwatch stopwatch = new Stopwatch();
            radius = radius / blockSize;

            string fname = Path.GetFileNameWithoutExtension(fpImgPath);

            Console.WriteLine("running LAS for {0}", fname);

            Image<Gray, Single> input = new Image<Gray, Single>(fpImgPath);

            Initialize(input);

            calGredient(ref input, ref Gx, ref Gy);

            calGxxGyyGxy(Gx, Gy, out Gxx, out Gyy, out Gxy);

            int windowsSize = (2 * (int)((2 * sigmaAvgWindows) + 0.5)) + 1;//20;
            Image<Gray, Single> avgWindow = new Image<Gray, Single>(windowsSize, windowsSize);
            calGaussianFilter(ref avgWindow, (float)sigmaAvgWindows);

            Convoluition(ref Gxx, ref Gxx, ref avgWindow);
            Convoluition(ref Gyy, ref Gyy, ref avgWindow);
            Convoluition(ref Gxy, ref Gxy, ref avgWindow);

            calPOF(Gxx, Gyy, Gxy, out POFx, out POFy);

            int GaussinWindowSize = (2 * (int)((2 * sigmaGaussian) + 0.5)) + 1;
            Image<Gray, Single> GaussianSmoothFilter = new Image<Gray, Single>(GaussinWindowSize, GaussinWindowSize);
            calGaussianFilter(ref GaussianSmoothFilter, (float)sigmaGaussian);

            Convoluition(ref POFx, ref POFx, ref GaussianSmoothFilter);
            Convoluition(ref POFy, ref POFy, ref GaussianSmoothFilter);

            Image<Gray, Single> POFxBlk = new Image<Gray, float>(1, 1);
            Image<Gray, Single> POFyBlk = new Image<Gray, float>(1, 1);
            pixelToBlock(POFx.Clone(), ref POFxBlk, blockSize);
            pixelToBlock(POFy.Clone(), ref POFyBlk, blockSize);

            Image<Gray, Single> inputBlk = input.Clone();
            pixelToBlock(input, ref inputBlk, blockSize);

            Image<Gray, Single> LASBlk_Value = new Image<Gray, float>(POFyBlk.Width, POFyBlk.Height, new Gray(0));
            Image<Gray, Single> LASBlk_Direction = new Image<Gray, float>(POFyBlk.Width, POFyBlk.Height, new Gray(0));

            for (int x = 0; x < POFxBlk.Width; x++)
            {
                for (int y = 0; y < POFxBlk.Height; y++)
                {
                    if (!float.IsNaN(POFxBlk.Data[y, x, 0]))
                    {
                        Point xy = new Point(x, y);
                        double direction = Double.NaN;
                        LASBlk_Value.Data[y, x, 0] = (float)calAOSforPointXY(POFxBlk, POFyBlk, xy, out direction);// calAOSforPointXYWithTable(POFxBlk, POFyBlk, xy, out direction, MapPixelForLASPath);
                        LASBlk_Direction.Data[y, x, 0] = (float)direction;
                    }
                    else
                    {
                        LASBlk_Value.Data[y, x, 0] = float.NaN;
                        LASBlk_Direction.Data[y, x, 0] = float.NaN;
                    }
                }
            }

            Image<Gray, Single> LASBlk_ValuePlot = LASBlk_Value.Clone();
            ScalingImage(ref LASBlk_ValuePlot);

            double[] LASMin, LASMax;
            Point[] LASPMin, LASPMax;
            LASBlk_Value.MinMax(out LASMin, out LASMax, out LASPMin, out LASPMax);

            using (var writer = new StreamWriter(Path.Combine(savePath, fname + "_MinMax.txt")))
            {
                writer.WriteLine("Min\tMax");
                writer.WriteLine("{0:R}\t{1:R}", LASMin[0], LASMax[0]);
            }

            LASBlk_ValuePlot.Save(savePath + fname + "_LASField.png");

            using (var writer = new StreamWriter(Path.Combine(savePath, fname + "_Direction.txt")))
            {
                for (int y = 0; y < LASBlk_Direction.Height; y++)
                {
                    for (int x = 0; x < LASBlk_Direction.Width; x++)
                    {
                        writer.Write("{0:R}\t", LASBlk_Direction.Data[y, x, 0]);
                    }
                    writer.Write("\n");
                }
            }

            Image<Gray, Single> LASBlk_Value60Max = LASBlk_Value.ThresholdBinaryInv(new Gray(LASMax[0] * 0.6), new Gray(1));
            Image<Gray, Single> LASBlk_Value60MaxPlot = LASBlk_Value60Max.Clone();
            ScalingImage(ref LASBlk_Value60MaxPlot, 255, 128);

            int m = 2;
            List<List<Point>> polygons = new List<List<Point>>();
            List<Point> centroids = new List<Point>();
            for (int x = 0; x < LASBlk_Value60Max.Width; x++)
            {
                for (int y = 0; y < LASBlk_Value60Max.Height; y++)
                {
                    if (LASBlk_Value60Max.Data[y, x, 0] == 1)
                    {
                        List<Point> polygonTmp = new List<Point>();
                        polygons.Add(polygonTmp);
                        Point centroidTmp = new Point();
                        floodFill(ref LASBlk_Value60Max, x, y, m, polygons[m - 2], ref centroidTmp);
                        ////discard area less  than 50 pixels
                        //if (polygons[m - 2].Count < 50)
                        //{
                        //    for (int c = 0; c < polygons[m - 2].Count; c++)
                        //    {
                        //        LASBlk_Value60Max.Data[polygons[m - 2][c].Y, polygons[m - 2][c].X, 0] = 0;
                        //    }
                        //    polygons.RemoveAt(m - 2);
                        //    m--;
                        //}
                        //else
                        //{
                        centroids.Add(centroidTmp);
                        //}
                        m++;
                    }
                }
            }

            LASBlk_Value60MaxPlot = LASBlk_Value60Max.Clone();
            ScalingImage(ref LASBlk_Value60MaxPlot, 255, 128);

            Image<Bgr, Single> LASBlk_Value60MaxBgrPlot = LASBlk_Value60MaxPlot.Convert<Bgr, Single>().Clone();
            Bgr centroidColor = new Bgr(Color.Red);
            for (int i = 0; i < centroids.Count; i++)
            {
                LASBlk_Value60MaxBgrPlot[centroids[i]] = centroidColor;
            }

            LASBlk_Value60MaxPlot.Save(savePath + fname + "_LASComp.png");
            LASBlk_Value60MaxBgrPlot.Save(savePath + fname + "_LASCompCentroid.png");

            input.Dispose();
            avgWindow.Dispose();
            GaussianSmoothFilter.Dispose();
            POFxBlk.Dispose();
            POFyBlk.Dispose();
            inputBlk.Dispose();
            LASBlk_Value.Dispose();
            LASBlk_Direction.Dispose();
            LASBlk_ValuePlot.Dispose();
            LASBlk_Value60Max.Dispose();
            LASBlk_Value60MaxPlot.Dispose();
            LASBlk_Value60MaxBgrPlot.Dispose();
            for (int i = 0; i < polygons.Count; i++)
            {
                polygons[i].Clear();
                polygons[i].TrimExcess();
                polygons[i] = null;
            }
            polygons.Clear();
            polygons.TrimExcess();
            polygons = null;
            centroids.Clear();
            centroids.TrimExcess();
            centroids = null;

            free();

            stopwatch.Stop();
            Console.WriteLine(fname + "'s LAS Field has been saved. Time elapsed: {0}", stopwatch.Elapsed);

        }

        private void calGredient(ref Image<Gray, Single> input, ref Image<Gray, Single> Gx, ref Image<Gray, Single> Gy)
        {
            Gx = input.Sobel(1, 0, 3);
            Gy = input.Sobel(0, 1, 3);
        }

        private void calGxxGyyGxy(Image<Gray, Single> Gx, Image<Gray, Single> Gy, out Image<Gray, Single> Gxx, out Image<Gray, Single> Gyy, out Image<Gray, Single> Gxy)
        {
            Image<Gray, Single> tmpGxx = new Image<Gray, float>(Gx.Width, Gy.Height);
            Image<Gray, Single> tmpGyy = new Image<Gray, float>(Gx.Width, Gy.Height);
            Image<Gray, Single> tmpGxy = new Image<Gray, float>(Gx.Width, Gy.Height);

            for (int x = 0; x < Gx.Width; x++)
            {
                for (int y = 0; y < Gx.Height; y++)
                {
                    tmpGxx.Data[y, x, 0] = Gx.Data[y, x, 0] * Gx.Data[y, x, 0]; //GxSign.Data[y, x, 0] * GxSign.Data[y, x, 0];
                    tmpGyy.Data[y, x, 0] = Gy.Data[y, x, 0] * Gy.Data[y, x, 0]; //GySign.Data[y, x, 0] * GySign.Data[y, x, 0];
                    tmpGxy.Data[y, x, 0] = Gx.Data[y, x, 0] * Gy.Data[y, x, 0]; //GxSign.Data[y, x, 0] * GySign.Data[y, x, 0];
                }
            }

            Gxx = tmpGxx.Clone();
            Gyy = tmpGyy.Clone();
            Gxy = tmpGxy.Clone();

            tmpGxx.Dispose();
            tmpGyy.Dispose();
            tmpGxy.Dispose();
        }

        private void calGaussianFilter(ref Image<Gray, Single> outFilter, float sigma, int derivativeMode = 0)
        {
            float sumTotal = 0;

            int kernelRadiusX = outFilter.Width / 2;
            int kernelRadiusY = outFilter.Height / 2;
            int EndkernelRadiusX = outFilter.Width / 2 + (outFilter.Width % 2 > 0 ? 1 : 0);
            int EndkernelRadiusY = outFilter.Height / 2 + (outFilter.Height % 2 > 0 ? 1 : 0);
            float distance;

            float calculatedEuler = (float)(1.0 / (2.0 * Math.PI * Math.Pow(sigma, 2)));

            for (int filterY = -kernelRadiusY;
                 filterY < EndkernelRadiusY; filterY++)
            {
                for (int filterX = -kernelRadiusX;
                    filterX < EndkernelRadiusX; filterX++)
                {
                    distance = (((filterX * filterX) + (filterY * filterY)) / (2 * (sigma * sigma)));

                    if (derivativeMode != 0)
                    {
                        if (derivativeMode == 1)
                        {
                            calculatedEuler = (float)(-filterX / (2.0 * Math.PI * Math.Pow(sigma, 4)));
                        }
                        else
                        {
                            calculatedEuler = (float)(-filterY / (2.0 * Math.PI * Math.Pow(sigma, 4)));
                        }
                    }


                    outFilter.Data[filterY + kernelRadiusY, filterX + kernelRadiusX, 0] = (float)(calculatedEuler * Math.Exp(-distance));


                    sumTotal += outFilter.Data[filterY + kernelRadiusY, filterX + kernelRadiusX, 0];
                }
            }

            for (int y = 0; y < outFilter.Height; y++)
            {
                for (int x = 0; x < outFilter.Width; x++)
                {
                    outFilter.Data[y, x, 0] = (float)(outFilter.Data[y, x, 0] * (1.0 / sumTotal));
                }
            }
        }

        private void Convoluition(ref Image<Gray, Single> input, ref Image<Gray, Single> output, ref Image<Gray, Single> mask)
        {
            Image<Gray, Single> tmpOutput = new Image<Gray, float>(input.Width, input.Height);
            tmpOutput.SetZero();

            for (int x = 0; x < input.Width; x++)
            {
                for (int y = 0; y < input.Height; y++)
                {
                    int StartBlockX = x - (mask.Width / 2);
                    int StartBlockY = y - (mask.Height / 2);
                    int EndBlockX = x + (mask.Width / 2) + (mask.Width % 2 > 0 ? 1 : 0);
                    int EndBlockY = y + (mask.Height / 2) + (mask.Height % 2 > 0 ? 1 : 0);

                    float countPix = 0;
                    float sumData = 0;
                    float avgData = 0;

                    for (int i = StartBlockX, iMask = 0; i < EndBlockX; i++, iMask++)
                    {
                        for (int j = StartBlockY, jMask = 0; j < EndBlockY; j++, jMask++)
                        {
                            if (i >= 0 && j >= 0 && i < input.Width && j < input.Height && !Double.IsNaN(input.Data[j, i, 0]))
                            {
                                sumData += input.Data[j, i, 0] * mask.Data[jMask, iMask, 0];
                                countPix++;
                            }
                        }
                    }

                    avgData = sumData / countPix;

                    tmpOutput.Data[y, x, 0] = avgData;
                }
            }

            output = tmpOutput.Clone();
            tmpOutput.Dispose();
        }

        private void calPOF(Image<Gray, Single> Gxx, Image<Gray, Single> Gyy, Image<Gray, Single> Gxy, out Image<Gray, Single> POFx, out Image<Gray, Single> POFy)
        {
            Image<Gray, Single> tmpPOFx = new Image<Gray, float>(Gxx.Width, Gxx.Height);
            Image<Gray, Single> tmpPOFy = new Image<Gray, float>(Gxx.Width, Gxx.Height);
            for (int x = 0; x < Gxx.Width; x++)
            {
                for (int y = 0; y < Gxx.Height; y++)
                {
                    tmpPOFx.Data[y, x, 0] = Gxx.Data[y, x, 0] - Gyy.Data[y, x, 0];
                    tmpPOFy.Data[y, x, 0] = 2 * Gxy.Data[y, x, 0];
                }
            }

            POFx = tmpPOFx.Clone();
            POFy = tmpPOFy.Clone();

            tmpPOFx.Dispose();
            tmpPOFy.Dispose();
        }

        private void pixelToBlock(Image<Gray, Single> pixelImg, ref Image<Gray, Single> outBlockImg, int blockSize)
        {
            int outputWidth = pixelImg.Width / blockSize;
            int outputHeight = pixelImg.Height / blockSize;
            if (pixelImg.Width % blockSize != 0)
            {
                outputWidth++;
            }
            if (pixelImg.Height % blockSize != 0)
            {
                outputHeight++;
            }

            Image<Gray, Single> blockImg = new Image<Gray, Single>(outputWidth, outputHeight, new Gray(0));

            for (int y = 0, j = 0; y < pixelImg.Height; y += blockSize, j++)
            {
                for (int x = 0, i = 0; x < pixelImg.Width; x += blockSize, i++)
                {
                    blockImg.Data[j, i, 0] = pixelImg.Data[y, x, 0];
                }
            }

            outBlockImg = blockImg.Clone();
            pixelImg.Dispose();
        }

        private double calAOSforPointXY(Image<Gray, Single> POFxBlk, Image<Gray, Single> POFyBlk, Point observed, out double direction)
        {
            double angularBandwidth = (2 * Math.PI) / numOfSectors;
            List<double> axialDirection = new List<double>();
            List<double> Sym_i = new List<double>();

            for (int i = 0; i < numOfSectors; i++)
            {
                axialDirection.Add(i * angularBandwidth);

                double sumCosAOS = 0;

                for (int j = 0; j < (numOfSectors / 2); j++)
                {
                    double startAngle_ipj = (i * angularBandwidth) + (j * angularBandwidth);
                    double endAngle_ipj = (i * angularBandwidth) + ((j + 1) * angularBandwidth);

                    double startAngle_imj_1 = (i * angularBandwidth) - ((j + 1) * angularBandwidth);
                    double endAngle_imj_1 = (i * angularBandwidth) - (j * angularBandwidth);

                    if (startAngle_ipj > (2 * Math.PI) + 0.0174533 || endAngle_ipj > (2 * Math.PI) + 0.0174533) //0.0174533 = 1 degree
                    {
                        startAngle_ipj = startAngle_ipj - (2 * Math.PI);
                        endAngle_ipj = endAngle_ipj - (2 * Math.PI);
                    }

                    if (startAngle_imj_1 < 0 || endAngle_imj_1 < 0)
                    {
                        startAngle_imj_1 = (2 * Math.PI) + startAngle_imj_1;
                        endAngle_imj_1 = (2 * Math.PI) + endAngle_imj_1;
                    }

                    double AOSx_ipj = 0;
                    double AOSy_ipj = 0;
                    double count_POF_ipj = 0;

                    double AOSx_imj_1 = 0;
                    double AOSy_imj_1 = 0;
                    double count_POF_imj_1 = 0;

                    for (int x = 0; x < POFxBlk.Width; x++)
                    {
                        for (int y = 0; y < POFxBlk.Height; y++)
                        {
                            if (!float.IsNaN(POFxBlk.Data[y, x, 0]))
                            {
                                Point xy = new Point(x, y);
                                double distance = distanceOf2Points(xy, observed);
                                double polarAngle = angleOf2Points(xy, observed);
                                polarAngle = polarAngle < 0 ? (2 * Math.PI) + polarAngle : polarAngle;

                                if (distance <= radius && polarAngle >= startAngle_ipj && polarAngle < endAngle_ipj)
                                {
                                    AOSx_ipj += POFxBlk.Data[y, x, 0];
                                    AOSy_ipj += POFyBlk.Data[y, x, 0];
                                    count_POF_ipj++;
                                }

                                if (distance <= radius && polarAngle >= startAngle_imj_1 && polarAngle < endAngle_imj_1)
                                {
                                    AOSx_imj_1 += POFxBlk.Data[y, x, 0];
                                    AOSy_imj_1 += POFyBlk.Data[y, x, 0];
                                    count_POF_imj_1++;
                                }
                            }
                        }
                    }

                    AOSx_ipj = AOSx_ipj / count_POF_ipj;
                    AOSy_ipj = AOSy_ipj / count_POF_ipj;

                    AOSx_imj_1 = AOSx_imj_1 / count_POF_imj_1;
                    AOSy_imj_1 = AOSy_imj_1 / count_POF_imj_1;

                    double AOS_ipj = Math.Atan2(AOSy_ipj, AOSx_ipj);
                    double AOS_imj_1 = Math.Atan2(AOSy_imj_1, AOSx_imj_1);

                    AOS_ipj = Double.IsNaN(AOS_ipj) ? 0 : AOS_ipj;
                    AOS_imj_1 = Double.IsNaN(AOS_imj_1) ? 0 : AOS_imj_1;

                    double cosAOS = Math.Cos(AOS_ipj + AOS_imj_1 - ((4 * i * 2 * Math.PI) / numOfSectors));

                    if (Double.IsNaN(cosAOS))
                    {
                        sumCosAOS += 0;
                    }
                    else
                    {
                        sumCosAOS += cosAOS;
                    }
                }

                double sym_i = (2 * sumCosAOS) / numOfSectors;
                Sym_i.Add(sym_i);
            }

            double LAS_xy = Sym_i.Max();
            int idxOfMaxValue = Sym_i.IndexOf(LAS_xy);
            direction = axialDirection[idxOfMaxValue];

            Sym_i.Clear();
            Sym_i.TrimExcess();
            Sym_i = null;

            return LAS_xy;
        }

        public void ScalingImage(ref Image<Gray, Single> inout, float maxValue = (float)255.0, float nanValue = (float)0.0)
        {
            double[] Min, Max;
            Point[] PMin, PMax;
            inout.MinMax(out Min, out Max, out PMin, out PMax);
            double divider = Max[0] - Min[0];
            for (int x = 0; x < inout.Width; x++)
            {
                for (int y = 0; y < inout.Height; y++)
                {
                    if (!Double.IsNaN(inout.Data[y, x, 0]))
                    {
                        Gray gValue = new Gray();
                        gValue.Intensity = (inout.Data[y, x, 0] - Min[0]) * (maxValue / divider);
                        inout[y, x] = gValue;
                    }
                    else
                    {
                        inout.Data[y, x, 0] = nanValue;
                    }
                }
            }
        }

        public static void floodFill(ref Image<Gray, float> ip, int x, int y, int label, List<Point> Polygon, ref Point centroid)
        {
            List<Point> s = new List<Point>(); // stack

            s.Add(new Point(x, y));
            int area = 0;
            int sum_x, sum_y, cx, cy;
            sum_x = sum_y = 0;

            while (s.Count > 0)
            {
                Point n = s.Last();
                s.RemoveAt(s.Count - 1);
                int u = n.X;
                int v = n.Y;

                if ((u >= 0) && (u < ip.Width) &&
                    (v >= 0) && (v < ip.Height) &&
                    ip.Data[v, u, 0] == 1)
                {
                    ip.Data[v, u, 0] = label;
                    Polygon.Add(n);
                    area++;
                    sum_x += u;
                    sum_y += v;

                    s.Add(new Point(u + 1, v));
                    s.Add(new Point(u, v + 1));
                    s.Add(new Point(u, v - 1));
                    s.Add(new Point(u - 1, v));
                }
            }
            cx = sum_x / area;
            cy = sum_y / area;
            centroid.X = cx;
            centroid.Y = cy;
        }

        public static double distanceOf2Points(PointF a, PointF b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public static double angleOf2Points(PointF a, PointF b)
        {
            return Math.Atan2(a.Y - b.Y, a.X - b.X);
        }
    }
}
