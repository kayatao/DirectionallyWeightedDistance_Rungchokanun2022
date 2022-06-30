using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//For PointF
using System.Drawing;
//For IO
using System.IO;
//For Image Utilities
using Emgu.CV;
using Emgu.CV.Structure;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class OF_Bazen2002
    {
        public Image<Gray, Single> Gx;
        public Image<Gray, Single> Gy;
        Image<Gray, Single> GxSign;
        Image<Gray, Single> GySign;
        public Image<Gray, Single> Gx2;
        public Image<Gray, Single> Gy2;
        public Image<Gray, Single> Gxy;
        Image<Gray, Single> GxxW;
        Image<Gray, Single> GyyW;
        Image<Gray, Single> GxyW;
        Image<Gray, Single> Phi;
        public Image<Gray, Single> Theta;
        public Image<Gray, Single> Coh;
        List<Point> NaNPoints;

        public OF_Bazen2002()
        {

        }

        ~OF_Bazen2002()
        {

        }

        public void free()
        {
            Gx.Dispose();
            Gy.Dispose();
            GxSign.Dispose();
            GySign.Dispose();
            Gx2.Dispose();
            Gy2.Dispose();
            Gxy.Dispose();
            GxxW.Dispose();
            GyyW.Dispose();
            GxyW.Dispose();
            Phi.Dispose();
            Theta.Dispose();
            Coh.Dispose();
            NaNPoints.Clear();
        }

        private void Initialize(Image<Gray, Single> input)
        {
            Gx = new Image<Gray, Single>(input.Width, input.Height);
            Gy = new Image<Gray, Single>(input.Width, input.Height);
            GxSign = new Image<Gray, Single>(input.Width, input.Height);
            GySign = new Image<Gray, Single>(input.Width, input.Height);
            Gx2 = new Image<Gray, Single>(input.Width, input.Height);
            Gy2 = new Image<Gray, Single>(input.Width, input.Height);
            Gxy = new Image<Gray, Single>(input.Width, input.Height);
            GxxW = new Image<Gray, Single>(input.Width, input.Height);
            GyyW = new Image<Gray, Single>(input.Width, input.Height);
            GxyW = new Image<Gray, Single>(input.Width, input.Height);
            Phi = new Image<Gray, Single>(input.Width, input.Height);
            Theta = new Image<Gray, Single>(input.Width, input.Height);
            Coh = new Image<Gray, Single>(input.Width, input.Height);
            NaNPoints = new List<Point>();
        }

        public void calOF(Image<Gray, Single> input, ref Image<Gray, Single> ThetaOutput, ref Image<Gray, Single> CohOutput, int windowsSize)
        {
            Initialize(input);
            calGredient(ref input, ref Gx, ref Gy, ref GxSign, ref GySign, false);

            for (int x = 0; x < GxSign.Width; x++)
            {
                for (int y = 0; y < GxSign.Height; y++)
                {
                    Gx2.Data[y, x, 0] = Gx.Data[y, x, 0] * Gx.Data[y, x, 0]; //GxSign.Data[y, x, 0] * GxSign.Data[y, x, 0];
                    Gy2.Data[y, x, 0] = Gy.Data[y, x, 0] * Gy.Data[y, x, 0]; //GySign.Data[y, x, 0] * GySign.Data[y, x, 0];
                    Gxy.Data[y, x, 0] = Gx.Data[y, x, 0] * Gy.Data[y, x, 0]; //GxSign.Data[y, x, 0] * GySign.Data[y, x, 0];
                }
            }

            Image<Gray, Single> OFWindow = new Image<Gray, Single>(windowsSize, windowsSize);
            calGaussianFilter(ref OFWindow, (float)5.0);

            Convoluition(ref Gx2, ref GxxW, ref OFWindow);
            Convoluition(ref Gy2, ref GyyW, ref OFWindow);
            Convoluition(ref Gxy, ref GxyW, ref OFWindow);

            for (int x = 0; x < Phi.Width; x++)
            {
                for (int y = 0; y < Phi.Height; y++)
                {
                    float angle = 0;
                    float Xterm = GxxW.Data[y, x, 0] - GyyW.Data[y, x, 0];
                    float Yterm = 2 * GxyW.Data[y, x, 0];

                    if (Xterm >= 0)
                    {
                        if (!Double.IsNaN(Math.Atan(Yterm / Xterm)))
                        {
                            angle = (float)Math.Atan(Yterm / Xterm);
                        }
                        else
                        {
                            angle = (float)Math.Atan(Yterm / Xterm);
                        }
                    }
                    else if (Xterm < 0 && Yterm >= 0)
                    {
                        if (!Double.IsNaN(Math.Atan(Yterm / Xterm)))
                        {
                            angle = (float)(Math.Atan(Yterm / Xterm) + Math.PI);
                        }
                        else
                        {
                            angle = 0 + (float)Math.PI;
                        }
                    }
                    else if (Xterm < 0 && Yterm < 0)
                    {
                        if (!Double.IsNaN(Math.Atan(Yterm / Xterm)))
                        {
                            angle = (float)(Math.Atan(Yterm / Xterm) - Math.PI);
                        }
                        else
                        {
                            angle = 0 - (float)Math.PI;
                        }
                    }

                    Phi.Data[y, x, 0] = (float)0.5 * angle;
                    if (Double.IsNaN(angle))
                    {
                        Point p = new Point(x, y);
                        NaNPoints.Add(p);
                    }

                    if (Phi.Data[y, x, 0] <= 0)
                    {
                        Theta.Data[y, x, 0] = Phi.Data[y, x, 0] + (float)(Math.PI / 2);
                    }
                    else if (Phi.Data[y, x, 0] > 0)
                    {
                        Theta.Data[y, x, 0] = Phi.Data[y, x, 0] - (float)(Math.PI / 2);
                    }
                    else
                    {
                        Theta.Data[y, x, 0] = Phi.Data[y, x, 0];
                    }

                    Coh.Data[y, x, 0] = (float)(Math.Sqrt(Math.Pow(GxxW.Data[y, x, 0] - GyyW.Data[y, x, 0], 2) + (4 * Math.Pow(GxyW.Data[y, x, 0], 2))) / (GxxW.Data[y, x, 0] + GyyW.Data[y, x, 0]));
                }
            }

            for (int i = 0; i < NaNPoints.Count; i++)
            {
                Gray NaNcolor = new Gray();
                NaNcolor.Intensity = Double.NaN;
                Theta[NaNPoints[i].Y, NaNPoints[i].X] = NaNcolor;
            }

            ThetaOutput = Theta.Clone();
            CohOutput = Coh.Clone();
        }

        private void calGredient(ref Image<Gray, Single> input, ref Image<Gray, Single> Gx, ref Image<Gray, Single> Gy, ref Image<Gray, Single> GxSign, ref Image<Gray, Single> GySign, bool bSmoothInput = false)
        {
            
            if (bSmoothInput)
            {
                Image<Gray, float> GaussianSmoothInput = new Image<Gray, float>(15, 15);
                calGaussianFilter(ref GaussianSmoothInput, (float)3.3);//3.3
                Matrix<float> mtrx = new Matrix<float>(GaussianSmoothInput.Width, GaussianSmoothInput.Height);
                GaussianSmoothInput.CopyTo(mtrx);
                ConvolutionKernelF kernel = new ConvolutionKernelF(mtrx, new Point(-1, -1));
                Image<Gray, Single> SmoothInput = input.Clone();
                SmoothInput = SmoothInput * kernel;

                Gx = SmoothInput.Sobel(1, 0, 3);
                Gy = SmoothInput.Sobel(0, 1, 3);

                GaussianSmoothInput.Dispose();
                mtrx.Dispose();
                kernel.Dispose();
                SmoothInput.Dispose();
            }
            else
            {
                Gx = input.Sobel(1, 0, 3);
                Gy = input.Sobel(0, 1, 3);
            }
            
            for (int x = 0; x < input.Width; x++)
            {
                for (int y = 0; y < input.Height; y++)
                {
                    float sign = 0;
                    Gray byCurrent = Gx[y, x];
                    if (byCurrent.Intensity < 0)
                    {
                        sign = -1;
                    }
                    else if (byCurrent.Intensity == 0)
                    {
                        sign = 0;
                    }
                    else if (byCurrent.Intensity > 0)
                    {
                        sign = 1;
                    }
                    GxSign.Data[y, x, 0] = sign * Gx.Data[y, x, 0];
                    GySign.Data[y, x, 0] = sign * Gy.Data[y, x, 0];
                }
            }
        }

        public void averageOF(ref Image<Gray, Single> theta, ref Image<Gray, Single> coh, ref Image<Gray, Single> outputTheta, ref Image<Gray, Single> outputCoh, int blocksize)
        {
            int outputWidth = theta.Width / blocksize;
            int outputHeight = theta.Height / blocksize;
            if (theta.Width % blocksize != 0)
            {
                outputWidth++;
            }
            if (theta.Height % blocksize != 0)
            {
                outputHeight++;
            }
            outputTheta = new Image<Gray, Single>(outputWidth, outputHeight);
            outputTheta.SetValue(0);
            outputCoh = new Image<Gray, Single>(outputWidth, outputHeight);
            outputCoh.SetValue(0);

            for (int x = blocksize / 2; x < theta.Width; x += blocksize)
            {
                for (int y = blocksize / 2; y < theta.Height; y += blocksize)
                {
                    int StartBlockX = x - (blocksize / 2);
                    int StartBlockY = y - (blocksize / 2);
                    int EndBlockX = StartBlockX + blocksize;
                    int EndBlockY = StartBlockY + blocksize;

                    float countPix = 0;
                    float sumThetaX = 0;
                    float sumThetaY = 0;
                    float sumCoh = 0;
                    float avgThetaX = 0;
                    float avgThetaY = 0;
                    float avgCoh = 0;
                    float avgTheta = 0;

                    for (int i = StartBlockX; i < EndBlockX; i++)
                    {
                        for (int j = StartBlockY; j < EndBlockY; j++)
                        {
                            if (i >= 0 && j >= 0 && i < theta.Width && j < theta.Height && !Double.IsNaN(theta.Data[j, i, 0]))
                            {
                                sumThetaX += (float)(coh.Data[j, i, 0] * Math.Cos(theta.Data[j, i, 0] * 2));
                                sumThetaY += (float)(coh.Data[j, i, 0] * Math.Sin(theta.Data[j, i, 0] * 2));

                                sumCoh += coh.Data[j, i, 0];
                                countPix++;
                            }
                        }
                    }

                    int outputX = StartBlockX / blocksize;
                    int outputY = StartBlockY / blocksize;
                    avgThetaX = sumThetaX / sumCoh; 
                    avgThetaY = sumThetaY / sumCoh; 
                    avgCoh = sumCoh / countPix;
                    avgTheta = (float)(Math.Atan2(avgThetaY, avgThetaX) / 2);

                    outputTheta.Data[outputY, outputX, 0] = avgTheta;
                    outputCoh.Data[outputY, outputX, 0] = avgCoh;
                }
            }
        }

        private void Convoluition(ref Image<Gray, Single> input, ref Image<Gray, Single> output, ref Image<Gray, Single> mask)
        {
            output.SetValue(0);

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

                    output.Data[y, x, 0] = avgData;
                }
            }
        }

        private void calGaussianFilter(ref Image<Gray, Single> outFilter, float sigma, int derivativeMode = 0)
        {
            float sumTotal = 0;


            int kernelRadiusX = outFilter.Width / 2;
            int kernelRadiusY = outFilter.Height / 2;
            int EndkernelRadiusX = outFilter.Width / 2 + (outFilter.Width % 2 > 0 ? 1 : 0);
            int EndkernelRadiusY = outFilter.Height / 2 + (outFilter.Height % 2 > 0 ? 1 : 0);
            float distance = 0;


            float calculatedEuler = (float)(1.0 / (2.0 * Math.PI * Math.Pow(sigma, 2)));


            for (int filterY = -kernelRadiusY;
                 filterY < EndkernelRadiusY; filterY++)
            {
                for (int filterX = -kernelRadiusX;
                    filterX < EndkernelRadiusX; filterX++)
                {
                    distance = (float)(((filterX * filterX) + (filterY * filterY)) / (2 * (sigma * sigma)));

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

    }
}
