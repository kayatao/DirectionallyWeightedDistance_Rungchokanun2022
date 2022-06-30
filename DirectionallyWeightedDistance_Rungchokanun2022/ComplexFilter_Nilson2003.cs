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
using Emgu.CV.CvEnum;
using Emgu.CV.UI;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class ComplexFilter_Nilson2003
    {
        public Image<Gray, Single> H1Re;
        public Image<Gray, Single> H1Im;
        public Image<Gray, Single> H2Re;
        public Image<Gray, Single> H2Im;
        public Image<Gray, Single> HGaussian;
        public float sigma = (float)1.5;

        public Image<Gray, Single> fx;
        public Image<Gray, Single> fy;

        public Image<Gray, Single>[] ZRek;
        public Image<Gray, Single>[] ZImk;

        public Image<Gray, Single>[] C1kRe;
        public Image<Gray, Single>[] C1kIm;
        public Image<Gray, Single>[] C2kRe;
        public Image<Gray, Single>[] C2kIm;

        public Image<Gray, Single>[] Mu1k;
        public Image<Gray, Single>[] Mu2k;
        public Image<Gray, Single>[] S1k;
        public Image<Gray, Single>[] S2k;

        public Image<Gray, Single>[] S1kDisp;
        public Image<Gray, Single>[] S2kDisp;

        public Image<Gray, Single>[] W1k;
        public Image<Gray, Single>[] W2k;

        public float ThetaCore;
        public float ThetaDelta;

        public Point[] /*core,*/ coreW;
        public Point[] /*delta,*/ deltaW;

        public bool Inited = false;

        double[][] S1kMin;
        double[][] S1kMax;
        Point[][] S1kPMin;
        Point[][] S1kPMax;
        double[][] S2kMin;
        double[][] S2kMax;
        Point[][] S2kPMin;
        Point[][] S2kPMax;

        double[][] S1kMinW;
        Point[][] S1kPMinW;
        double[][] S2kMinW;
        Point[][] S2kPMinW;

        public double[][] S1kMaxW;
        public Point[][] S1kPMaxW;
        public double[][] S2kMaxW;
        public Point[][] S2kPMaxW;

        Image<Gray, Single> Gaussian;
        Image<Gray, Single> GaussianX;
        Image<Gray, Single> GaussianY;
        Image<Gray, Single> GaussianSmoothZ3;

        Image<Gray, Single> Z3Mag;
        Image<Gray, Single> Z3MagG1;

        Image<Gray, Single> C13Z3G1Re;
        Image<Gray, Single> C13Z3G1Im;
        Image<Gray, Single> C13G2Re;
        Image<Gray, Single> C13G2Im;

        Image<Gray, Single> C23Z3G1Re;
        Image<Gray, Single> C23Z3G1Im;
        Image<Gray, Single> C23G2Re;
        Image<Gray, Single> C23G2Im;

        Image<Gray, Single> G2;

        Image<Gray, Single> C13WeightedRe;
        Image<Gray, Single> C13WeightedIm;
        Image<Gray, Single> C23WeightedRe;
        Image<Gray, Single> C23WeightedIm;

        public ComplexFilter_Nilson2003()
        {

        }

        ~ComplexFilter_Nilson2003()
        {

        }

        public void Initialize()
        {
            Inited = true;
            H1Re = new Image<Gray, Single>(9, 9);
            H1Im = new Image<Gray, Single>(9, 9);
            H2Re = new Image<Gray, Single>(9, 9);
            H2Im = new Image<Gray, Single>(9, 9);
            HGaussian = new Image<Gray, Single>(9, 9);
            calGaussianFilter(ref HGaussian, (float)sigma);

            genComplexFilter(ref H1Re, ref H1Im, (float)1.0, ref HGaussian);
            genComplexFilter(ref H2Re, ref H2Im, (float)-1.0, ref HGaussian);

            Gaussian = new Image<Gray, Single>(5, 5);
            GaussianX = new Image<Gray, Single>(5, 5);
            GaussianY = new Image<Gray, Single>(5, 5);
            GaussianSmoothZ3 = new Image<Gray, Single>(5, 5);

            calGaussianFilter(ref Gaussian, (float)0.8);
            calGaussianFilter2D(ref GaussianSmoothZ3, (float)1.5, (float)1.5);
            calGaussianFilter(ref GaussianX, (float)0.8, 1);
            calGaussianFilter(ref GaussianY, (float)0.8, 2);

            ZRek = new Image<Gray, float>[4];
            ZImk = new Image<Gray, float>[4];

            C1kRe = new Image<Gray, float>[4];
            C1kIm = new Image<Gray, float>[4];
            C2kRe = new Image<Gray, float>[4];
            C2kIm = new Image<Gray, float>[4];

            Mu1k = new Image<Gray, float>[4];
            Mu2k = new Image<Gray, float>[4];
            S1k = new Image<Gray, float>[4];
            S2k = new Image<Gray, float>[4];

            W1k = new Image<Gray, float>[4];
            W2k = new Image<Gray, float>[4];

            S1kDisp = new Image<Gray, float>[4];
            S2kDisp = new Image<Gray, float>[4];

            S1kMin = new double[4][];
            S1kMax = new double[4][];
            S1kPMin = new Point[4][];
            S1kPMax = new Point[4][];
            S2kMin = new double[4][];
            S2kMax = new double[4][];
            S2kPMin = new Point[4][];
            S2kPMax = new Point[4][];

            S1kMinW = new double[4][];
            S1kMaxW = new double[4][];
            S1kPMinW = new Point[4][];
            S1kPMaxW = new Point[4][];
            S2kMinW = new double[4][];
            S2kMaxW = new double[4][];
            S2kPMinW = new Point[4][];
            S2kPMaxW = new Point[4][];
        }

        public void free()
        {
            if (Inited)
            {
                Inited = false;
                H1Re.Dispose();
                H1Im.Dispose();
                H2Re.Dispose();
                H2Im.Dispose();
                HGaussian.Dispose();

                fx.Dispose();
                fy.Dispose();

                for (int i = 0; i < ZRek.Length; i++)
                {
                    ZRek[i].Dispose();
                    ZImk[i].Dispose();
                    C1kRe[i].Dispose();
                    C1kIm[i].Dispose();
                    C2kRe[i].Dispose();
                    C2kIm[i].Dispose();
                    Mu1k[i].Dispose();
                    Mu2k[i].Dispose();
                    S1k[i].Dispose();
                    S2k[i].Dispose();
                    S1kDisp[i].Dispose();
                    S2kDisp[i].Dispose();
                    W1k[i].Dispose();
                    W2k[i].Dispose();
                }

                Gaussian.Dispose();
                GaussianX.Dispose();
                GaussianY.Dispose();
                GaussianSmoothZ3.Dispose();

                Z3Mag.Dispose();
                Z3MagG1.Dispose();

                C13Z3G1Re.Dispose();
                C13Z3G1Im.Dispose();
                C13G2Re.Dispose();
                C13G2Im.Dispose();

                C23Z3G1Re.Dispose();
                C23Z3G1Im.Dispose();
                C23G2Re.Dispose();
                C23G2Im.Dispose();

                G2.Dispose();

                C13WeightedRe.Dispose();
                C13WeightedIm.Dispose();
                C23WeightedRe.Dispose();
                C23WeightedIm.Dispose();
            }
        }

        public void calGaussianFilter(ref Image<Gray, Single> outFilter, float sigma, int derivativeMode = 0)
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
                    distance = (float)((((float)filterX * (float)filterX) + ((float)filterY * (float)filterY)) / ((float)2.0 * (sigma * sigma)));

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

        public void calGaussianFilter2D(ref Image<Gray, Single> outFilter, float sigmaX, float sigmaY)
        {
            float sumTotal = 0;

            float MuX = outFilter.Width / 2;
            float MuY = outFilter.Height / 2;

            float calculatedEuler = (float)(1.0 / (2.0 * Math.PI * sigmaX * sigmaY));

            for (int filterY = 0; filterY < outFilter.Height; filterY++)
            {
                for (int filterX = 0; filterX < outFilter.Width; filterX++)
                {
                    float distance = (float)(((Math.Pow(filterX - MuX, 2) / (2 * Math.Pow(sigmaX, 2))) + (Math.Pow(filterY - MuY, 2) / (2 * Math.Pow(sigmaY, 2)))));

                    outFilter.Data[filterY, filterX, 0] = (float)(calculatedEuler * Math.Exp(-distance));

                    sumTotal += outFilter.Data[filterY, filterX, 0];
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

        private void genComplexFilter(ref Image<Gray, Single> hRe, ref Image<Gray, Single> hIm, float order, ref Image<Gray, Single> Gaussian/*, float sigmaOfGaussian*/)
        {
            int kernelRadiusX = Gaussian.Width / 2;
            int kernelRadiusY = Gaussian.Height / 2;
            int EndkernelRadiusX = Gaussian.Width / 2 + (Gaussian.Width % 2 > 0 ? 1 : 0);
            int EndkernelRadiusY = Gaussian.Height / 2 + (Gaussian.Height % 2 > 0 ? 1 : 0);

            for (int x = -kernelRadiusX; x < EndkernelRadiusX; x++)
            {
                for (int y = -kernelRadiusY; y < EndkernelRadiusY; y++)
                {
                    //include Magnitude
                    //Complex h1 = new Complex(x, y);
                    //h1 = Complex.Pow(h1, order);
                    //hRe.Data[y + kernelRadiusY, x + kernelRadiusX, 0] = (float)h1.Real * Gaussian.Data[y + kernelRadiusY, x + kernelRadiusX, 0];
                    //hIm.Data[y + kernelRadiusY, x + kernelRadiusX, 0] = (float)h1.Imaginary * Gaussian.Data[y + kernelRadiusY, x + kernelRadiusX, 0];

                    //ignore magnitude
                    float tmpHRe = x;
                    float tmpHIm = y;
                    if (order == -1)
                    {
                        tmpHIm = -tmpHIm;
                    }
                    hRe.Data[y + kernelRadiusY, x + kernelRadiusX, 0] = (float)tmpHRe * Gaussian.Data[y + kernelRadiusY, x + kernelRadiusX, 0];
                    hIm.Data[y + kernelRadiusY, x + kernelRadiusX, 0] = (float)tmpHIm * Gaussian.Data[y + kernelRadiusY, x + kernelRadiusX, 0];
                }
            }
        }

        public void calFilterResponse(Image<Gray, Single> input, int windowsSize)
        {
            #region Calculate tensor structure and create pyramid 4 levels

            for (int i = 0; i < ZRek.Length; i++)
            {
                if (i == 0)
                {
                    fx = new Image<Gray, float>(input.Width, input.Height);
                    fy = new Image<Gray, float>(input.Width, input.Height);

                    ConvolutionIgnoreOutOfMask(ref input, ref fx, ref GaussianX);
                    ConvolutionIgnoreOutOfMask(ref input, ref fy, ref GaussianY);

                    ZRek[i] = new Image<Gray, Single>(fx.Width, fx.Height, new Gray(0));
                    ZImk[i] = new Image<Gray, Single>(fx.Width, fx.Height, new Gray(0));
                    for (int x = 0; x < fx.Width; x++)
                    {
                        for (int y = 0; y < fx.Height; y++)
                        {

                            ZRek[i].Data[y, x, 0] = 0;
                            ZImk[i].Data[y, x, 0] = 0;

                            int StartBlockX = x - (windowsSize / 2);
                            int StartBlockY = y - (windowsSize / 2);
                            int EndBlockX = x + (windowsSize / 2) + (windowsSize % 2 > 0 ? 1 : 0);
                            int EndBlockY = y + (windowsSize / 2) + (windowsSize % 2 > 0 ? 1 : 0);

                            int countPix = 0;

                            for (int u = StartBlockX, iMask = 0; u < EndBlockX; u++, iMask++)
                            {
                                for (int v = StartBlockY, jMask = 0; v < EndBlockY; v++, jMask++)
                                {
                                    int xx = 0, yy = 0;
                                    MirrorCoordinate(u, v, ref xx, ref yy, ZRek[i].Width, ZRek[i].Height);
                                    float Gxx = (float)Math.Pow(fx.Data[yy, xx, 0], 2);
                                    float Gyy = (float)Math.Pow(fy.Data[yy, xx, 0], 2);
                                    float Gxy2 = 2 * fx.Data[yy, xx, 0] * fy.Data[yy, xx, 0];
                                    float ZRe = Gxx - Gyy;
                                    float ZIm = Gxy2;
                                    ZRek[i].Data[y, x, 0] += ZRe;
                                    ZImk[i].Data[y, x, 0] += ZIm;
                                    countPix++;
                                    //ZRek[i].Data[y, x, 0] += (float)(Math.Pow(fx.Data[yy, xx, 0], 2) - Math.Pow(fy.Data[yy, xx, 0], 2));
                                    //ZImk[i].Data[y, x, 0] += (float)(2.0 * fx.Data[yy, xx, 0] * fy.Data[yy, xx, 0]);
                                }
                            }
                            ZRek[i].Data[y, x, 0] = ZRek[i].Data[y, x, 0] / countPix;
                            ZImk[i].Data[y, x, 0] = ZImk[i].Data[y, x, 0] / countPix;
                        }
                    }
                }
                else
                {
                    ZRek[i] = ZRek[i - 1].Clone();
                    ZImk[i] = ZImk[i - 1].Clone();
                    ConvolutionIgnoreOutOfMask(ref ZRek[i - 1], ref ZRek[i], ref Gaussian);
                    ConvolutionIgnoreOutOfMask(ref ZImk[i - 1], ref ZImk[i], ref Gaussian);
                    Image<Gray, float> tmpX = new Image<Gray, float>((int)(ZRek[i].Width / 2.0 + 0.5), (int)(ZRek[i].Height / 2.0 + 0.5));
                    Image<Gray, float> tmpY = new Image<Gray, float>((int)(ZImk[i].Width / 2.0 + 0.5), (int)(ZImk[i].Height / 2.0 + 0.5));
                    CvInvoke.PyrDown(ZRek[i], tmpX);
                    CvInvoke.PyrDown(ZImk[i], tmpY);
                    ZRek[i] = tmpX.Clone();
                    ZImk[i] = tmpY.Clone();
                    tmpX.Dispose();
                    tmpY.Dispose();
                }

            }

            #endregion

            for (int i = 0; i < C1kRe.Length; i++)
            {
                C1kRe[i] = new Image<Gray, float>(ZRek[i].Width, ZRek[i].Height);
                C1kIm[i] = new Image<Gray, float>(ZImk[i].Width, ZImk[i].Height);
                C2kRe[i] = new Image<Gray, float>(ZRek[i].Width, ZRek[i].Height);
                C2kIm[i] = new Image<Gray, float>(ZImk[i].Width, ZImk[i].Height);
                ComplexConvoluition(ref ZRek[i], ref ZImk[i], ref C1kRe[i], ref C1kIm[i], ref H1Re, ref H1Im);
                ComplexConvoluition(ref ZRek[i], ref ZImk[i], ref C2kRe[i], ref C2kIm[i], ref H2Re, ref H2Im);
            }

            #region Adjust weight of filter response in level 3

            C13Z3G1Re = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C13Z3G1Im = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C13G2Re = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C13G2Im = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);

            C23Z3G1Re = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C23Z3G1Im = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C23G2Re = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C23G2Im = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);

            G2 = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);

            C13WeightedRe = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C13WeightedIm = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C23WeightedRe = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);
            C23WeightedIm = new Image<Gray, float>(C1kRe[3].Width, C1kRe[3].Height);

            calGaussianFilter2D(ref G2, (float)7.0, (float)11.7);

            Z3Mag = new Image<Gray, float>(ZRek[3].Width, ZRek[3].Height);
            Z3MagG1 = new Image<Gray, float>(ZRek[3].Width, ZRek[3].Height);

            for (int x = 0; x < Z3Mag.Width; x++)
            {
                for (int y = 0; y < Z3Mag.Height; y++)
                {
                    Z3Mag.Data[y, x, 0] = (float)Math.Sqrt((ZRek[3].Data[y, x, 0] * ZRek[3].Data[y, x, 0]) + (ZImk[3].Data[y, x, 0] * ZImk[3].Data[y, x, 0]));
                }
            }

            ScalingImage(ref G2, (float)1.0);

            Convolution(ref Z3Mag, ref Z3MagG1, ref GaussianSmoothZ3);
            ScalingImage(ref Z3MagG1, (float)1.0);

            for (int x = 0; x < Z3Mag.Width; x++)
            {
                for (int y = 0; y < Z3Mag.Height; y++)
                {
                    C13Z3G1Re.Data[y, x, 0] = Z3MagG1.Data[y, x, 0] * C1kRe[3].Data[y, x, 0];
                    C13Z3G1Im.Data[y, x, 0] = Z3MagG1.Data[y, x, 0] * C1kIm[3].Data[y, x, 0];
                    C13G2Re.Data[y, x, 0] = G2.Data[y, x, 0] * C1kRe[3].Data[y, x, 0];
                    C13G2Im.Data[y, x, 0] = G2.Data[y, x, 0] * C1kIm[3].Data[y, x, 0];

                    C23Z3G1Re.Data[y, x, 0] = Z3MagG1.Data[y, x, 0] * C2kRe[3].Data[y, x, 0];
                    C23Z3G1Im.Data[y, x, 0] = Z3MagG1.Data[y, x, 0] * C2kIm[3].Data[y, x, 0];
                    C23G2Re.Data[y, x, 0] = G2.Data[y, x, 0] * C2kRe[3].Data[y, x, 0];
                    C23G2Im.Data[y, x, 0] = G2.Data[y, x, 0] * C2kIm[3].Data[y, x, 0];

                    C1kRe[3].Data[y, x, 0] = (C13Z3G1Re.Data[y, x, 0] + C13G2Re.Data[y, x, 0]) / 2;
                    C1kIm[3].Data[y, x, 0] = (C13Z3G1Im.Data[y, x, 0] + C13G2Im.Data[y, x, 0]) / 2;

                    C2kRe[3].Data[y, x, 0] = (C23Z3G1Re.Data[y, x, 0] + C23G2Re.Data[y, x, 0]) / 2;
                    C2kIm[3].Data[y, x, 0] = (C23Z3G1Im.Data[y, x, 0] + C23G2Im.Data[y, x, 0]) / 2;
                }
            }

            #endregion

            #region Calculate Magnitude response S1k and S2k

            for (int i = 0; i < Mu1k.Length; i++)
            {
                Mu1k[i] = new Image<Gray, Single>(ZRek[i].Width, ZRek[i].Height);
                Mu2k[i] = new Image<Gray, Single>(ZRek[i].Width, ZRek[i].Height);
                for (int x = 0; x < ZRek[i].Width; x++)
                {
                    for (int y = 0; y < ZRek[i].Height; y++)
                    {
                        Mu1k[i].Data[y, x, 0] = (float)Math.Sqrt((C1kRe[i].Data[y, x, 0] * C1kRe[i].Data[y, x, 0]) + (C1kIm[i].Data[y, x, 0] * C1kIm[i].Data[y, x, 0]));
                        Mu2k[i].Data[y, x, 0] = (float)Math.Sqrt((C2kRe[i].Data[y, x, 0] * C2kRe[i].Data[y, x, 0]) + (C2kIm[i].Data[y, x, 0] * C2kIm[i].Data[y, x, 0]));
                    }
                }

                ScalingImage(ref Mu1k[i], (float)1.0);
                ScalingImage(ref Mu2k[i], (float)1.0);

                S1k[i] = new Image<Gray, Single>(ZRek[i].Width, ZRek[i].Height);
                S2k[i] = new Image<Gray, Single>(ZRek[i].Width, ZRek[i].Height);
                for (int x = 0; x < ZRek[i].Width; x++)
                {
                    for (int y = 0; y < ZRek[i].Height; y++)
                    {
                        S1k[i].Data[y, x, 0] = Mu1k[i].Data[y, x, 0] * (1 - Mu2k[i].Data[y, x, 0]);
                        S2k[i].Data[y, x, 0] = Mu2k[i].Data[y, x, 0] * (1 - Mu1k[i].Data[y, x, 0]);
                    }
                }
            }

            #endregion

            #region Searching maximum response

            for (int i = 0; i < S1k.Length; i++)
            {
                S1k[i].MinMax(out S1kMin[i], out S1kMax[i], out S1kPMin[i], out S1kPMax[i]);
                S2k[i].MinMax(out S2kMin[i], out S2kMax[i], out S2kPMin[i], out S2kPMax[i]);

                S1kMinW[i] = S1kMin[i];
                S1kMaxW[i] = S1kMax[i];
                S1kPMinW[i] = S1kPMin[i];
                S1kPMaxW[i] = S1kPMax[i];

                S2kMinW[i] = S2kMin[i];
                S2kMaxW[i] = S2kMax[i];
                S2kPMinW[i] = S2kPMin[i];
                S2kPMaxW[i] = S2kPMax[i];

                S1kDisp[i] = S1k[i].Clone();
                S2kDisp[i] = S2k[i].Clone();
                ScalingImage0To1(ref S1kDisp[i]);
                ScalingImage0To1(ref S2kDisp[i]);
            }

            for (int i = S1k.Length - 1; i > 0; i--)
            {
                int multiplier = (int)Math.Pow(2, 1);
                int offset = 2;
                int x1ktox1k_1 = S1kPMax[i][0].X * multiplier + offset;
                int y1ktoy1k_1 = S1kPMax[i][0].Y * multiplier + offset;

                int x2ktox2k_1 = S2kPMax[i][0].X * multiplier + offset;
                int y2ktoy2k_1 = S2kPMax[i][0].Y * multiplier + offset;

                float maxS1k_1 = -9999, maxS2k_1 = -9999;
                int x1k_1 = 0, y1k_1 = 0, x2k_1 = 0, y2k_1 = 0;

                for (int x = x1ktox1k_1 - 6; x <= x1ktox1k_1 + 6; x++)
                {
                    for (int y = y1ktoy1k_1 - 6; y <= y1ktoy1k_1 + 6; y++)
                    {
                        if (x >= 0 && x < S1k[i - 1].Width && y >= 0 && y < S1k[i - 1].Height)
                        {
                            if (maxS1k_1 < S1k[i - 1].Data[y, x, 0])
                            {
                                maxS1k_1 = S1k[i - 1].Data[y, x, 0];
                                x1k_1 = x;
                                y1k_1 = y;
                            }
                        }
                    }
                }

                for (int x = x2ktox2k_1 - 6; x <= x2ktox2k_1 + 6; x++)
                {
                    for (int y = y2ktoy2k_1 - 6; y <= y2ktoy2k_1 + 6; y++)
                    {
                        if (x >= 0 && x < S2k[i - 1].Width && y >= 0 && y < S2k[i - 1].Height)
                        {
                            if (maxS2k_1 < S2k[i - 1].Data[y, x, 0])
                            {
                                maxS2k_1 = S2k[i - 1].Data[y, x, 0];
                                x2k_1 = x;
                                y2k_1 = y;
                            }
                        }
                    }
                }

                S1kMaxW[i - 1][0] = maxS1k_1;
                S1kPMaxW[i - 1][0].X = x1k_1;
                S1kPMaxW[i - 1][0].Y = y1k_1;

                S2kMaxW[i - 1][0] = maxS2k_1;
                S2kPMaxW[i - 1][0].X = x2k_1;
                S2kPMaxW[i - 1][0].Y = y2k_1;

            }

            #endregion

            #region Calculate orientation of singular points

            for (int k = 0; k < C1kRe.Length; k++)
            {
                W1k[k] = new Image<Gray, float>(C1kRe[k].Width, C1kRe[k].Height);
                W2k[k] = new Image<Gray, float>(C1kRe[k].Width, C1kRe[k].Height);

                int weightMaskSize = 3;

                int StartBlockX = S1kPMaxW[k][0].X - (weightMaskSize / 2);
                int StartBlockY = S1kPMaxW[k][0].Y - (weightMaskSize / 2);
                int EndBlockX = S1kPMaxW[k][0].X + (weightMaskSize / 2) + (weightMaskSize % 2 > 0 ? 1 : 0);
                int EndBlockY = S1kPMaxW[k][0].Y + (weightMaskSize / 2) + (weightMaskSize % 2 > 0 ? 1 : 0);
                float denominator1 = 0;
                float denominator2 = 0;

                for (int i = StartBlockX, iMask = 0; i < EndBlockX; i++, iMask++)
                {
                    for (int j = StartBlockY, jMask = 0; j < EndBlockY; j++, jMask++)
                    {
                        if (i >= 0 && j >= 0 && i < C1kRe[k].Width && j < C1kRe[k].Height)
                        {
                            float magnitude1 = (float)(Math.Sqrt(Math.Pow(C1kRe[k].Data[j, i, 0], 2) + Math.Pow(C1kIm[k].Data[j, i, 0], 2)));
                            denominator1 += magnitude1;
                        }
                    }
                }

                StartBlockX = S2kPMaxW[k][0].X - (weightMaskSize / 2);
                StartBlockY = S2kPMaxW[k][0].Y - (weightMaskSize / 2);
                EndBlockX = S2kPMaxW[k][0].X + (weightMaskSize / 2) + (weightMaskSize % 2 > 0 ? 1 : 0);
                EndBlockY = S2kPMaxW[k][0].Y + (weightMaskSize / 2) + (weightMaskSize % 2 > 0 ? 1 : 0);
                for (int i = StartBlockX, iMask = 0; i < EndBlockX; i++, iMask++)
                {
                    for (int j = StartBlockY, jMask = 0; j < EndBlockY; j++, jMask++)
                    {
                        if (i >= 0 && j >= 0 && i < C2kRe[k].Width && j < C2kRe[k].Height)
                        {
                            float magnitude2 = (float)(Math.Sqrt(Math.Pow(C2kRe[k].Data[j, i, 0], 2) + Math.Pow(C2kIm[k].Data[j, i, 0], 2)));
                            denominator2 += magnitude2;
                        }
                    }
                }

                for (int x = 0; x < C1kRe[k].Width; x++)
                {
                    for (int y = 0; y < C1kRe[k].Height; y++)
                    {

                        float numerator1 = (float)(Math.Sqrt(Math.Pow(C1kRe[k].Data[y, x, 0], 2) + Math.Pow(C1kIm[k].Data[y, x, 0], 2)));
                        float avgData1 = 0;

                        float numerator2 = (float)(Math.Sqrt(Math.Pow(C2kRe[k].Data[y, x, 0], 2) + Math.Pow(C2kIm[k].Data[y, x, 0], 2)));
                        float avgData2 = 0;

                        if (denominator1 != 0)
                        {
                            avgData1 = numerator1 / denominator1;
                        }

                        if (denominator2 != 0)
                        {
                            avgData2 = numerator2 / denominator2;
                        }

                        W1k[k].Data[y, x, 0] = avgData1;
                        W2k[k].Data[y, x, 0] = avgData2;
                    }
                }
                
            }

            float[] W1kC1kRe = new float[4];
            float[] W1kC1kIm = new float[4];
            float[] W2kC2kRe = new float[4];
            float[] W2kC2kIm = new float[4];

            int neighbourhoodSize = 3;

            for (int k = 0; k < C1kRe.Length; k++)
            {
                int max1PX = S1kPMaxW[k][0].X;
                int max1PY = S1kPMaxW[k][0].Y;

                int max2PX = S2kPMaxW[k][0].X;
                int max2PY = S2kPMaxW[k][0].Y;

                W1kC1kRe[k] = 0;
                W1kC1kIm[k] = 0;
                W2kC2kRe[k] = 0;
                W2kC2kIm[k] = 0;

                int countpix = 0;
                for (int x = max1PX - (neighbourhoodSize / 2); x <= max1PX + (neighbourhoodSize / 2); x++)
                {
                    for (int y = max1PY - (neighbourhoodSize / 2); y <= max1PY + (neighbourhoodSize / 2); y++)
                    {
                        if (x >= 0 && y >= 0 && x < W1k[k].Width && y < W1k[k].Height)
                        {
                            W1kC1kRe[k] += W1k[k].Data[y, x, 0] * C1kRe[k].Data[y, x, 0];
                            W1kC1kIm[k] += W1k[k].Data[y, x, 0] * C1kIm[k].Data[y, x, 0];
                            countpix++;
                        }
                    }
                }

                countpix = 0;
                for (int x = max2PX - (neighbourhoodSize / 2); x <= max2PX + (neighbourhoodSize / 2); x++)
                {
                    for (int y = max2PY - (neighbourhoodSize / 2); y <= max2PY + (neighbourhoodSize / 2); y++)
                    {
                        if (x >= 0 && y >= 0 && x < W1k[k].Width && y < W1k[k].Height)
                        {
                            W2kC2kRe[k] += W2k[k].Data[y, x, 0] * C2kRe[k].Data[y, x, 0];
                            W2kC2kIm[k] += W2k[k].Data[y, x, 0] * C2kIm[k].Data[y, x, 0];
                            countpix++;
                        }
                    }
                }
            }

            float sumRealTerm1 = W1kC1kRe[2] + W1kC1kRe[3];
            float sumImTerm1 = W1kC1kIm[2] + W1kC1kIm[3];

            float sumRealTerm2 = W2kC2kRe[2] + W2kC2kRe[3];
            float sumImTerm2 = W2kC2kIm[2] + W2kC2kIm[3];

            ThetaCore = (float)(Math.Atan2(sumImTerm1, sumRealTerm1));
            ThetaDelta = (float)((Math.Atan2(sumImTerm2, sumRealTerm2) / 3.0));


            #endregion

        }

        public void ConvolutionIgnoreOutOfMask(ref Image<Gray, Single> input, ref Image<Gray, Single> output, ref Image<Gray, Single> mask)
        {
            Image<Gray, Single> tmpOut = new Image<Gray, float>(input.Width, input.Height);

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

                    if (StartBlockX >= 0 && StartBlockY >= 0 && EndBlockX <= input.Width && EndBlockY <= input.Height)
                    {
                        for (int i = StartBlockX, iMask = 0; i < EndBlockX; i++, iMask++)
                        {
                            for (int j = StartBlockY, jMask = 0; j < EndBlockY; j++, jMask++)
                            {
                                sumData += input.Data[j, i, 0] * mask.Data[jMask, iMask, 0];
                                countPix++;
                            }
                        }

                        avgData = sumData / countPix;

                    }

                    tmpOut.Data[y, x, 0] = avgData;
                }
            }
            int strROIX = 0 + (mask.Width / 2);
            int strROIY = 0 + (mask.Height / 2);
            int widthROI = input.Width - (mask.Width / 2) - strROIX;
            int heightROI = input.Height - (mask.Height / 2) - strROIY;
            tmpOut.ROI = new Rectangle(strROIX, strROIY, widthROI, heightROI);
            tmpOut = tmpOut.Copy();

            output = tmpOut.Clone();
            tmpOut.Dispose();
        }

        public void MirrorCoordinate(int x, int y, ref int new_x, ref int new_y, int width, int height)
        {
            if (x >= 0 && x < width) new_x = x;
            else if (x < 0) new_x = Math.Abs(x);
            else new_x = x - width;

            if (y >= 0 && y < height) new_y = y;
            else if (y < 0) new_y = Math.Abs(y);
            else new_y = y - height;
        }

        public void ComplexConvoluition(ref Image<Gray, Single> iRe, ref Image<Gray, Single> iIm, ref Image<Gray, Single> oRe, ref Image<Gray, Single> oIm, ref Image<Gray, Single> filterRe, ref Image<Gray, Single> filterIm)
        {
            Image<Gray, Single> tmpOutRe = new Image<Gray, float>(iRe.Width, iRe.Height);
            Image<Gray, Single> tmpOutIm = new Image<Gray, float>(iIm.Width, iIm.Height);

            for (int x = 0; x < iRe.Width; x++)
            {
                for (int y = 0; y < iRe.Height; y++)
                {
                    int StartBlockX = x - (filterRe.Width / 2);
                    int StartBlockY = y - (filterRe.Height / 2);
                    int EndBlockX = x + (filterRe.Width / 2) + (filterRe.Width % 2 > 0 ? 1 : 0);
                    int EndBlockY = y + (filterRe.Height / 2) + (filterRe.Height % 2 > 0 ? 1 : 0);

                    float sumDataRe = 0;
                    float sumDataIm = 0;

                    for (int i = StartBlockX, iMask = 0; i < EndBlockX; i++, iMask++)
                    {
                        for (int j = StartBlockY, jMask = 0; j < EndBlockY; j++, jMask++)
                        {
                            if (i >= 0 && j >= 0 && i < iRe.Width && j < iRe.Height)
                            {
                                float ZRe = iRe.Data[j, i, 0];
                                float ZIm = iIm.Data[j, i, 0];
                                sumDataRe += (ZRe * filterRe.Data[jMask, iMask, 0]) - (ZIm * filterIm.Data[jMask, iMask, 0]);
                                sumDataIm += (ZIm * filterRe.Data[jMask, iMask, 0]) + (ZRe * filterIm.Data[jMask, iMask, 0]);
                            }
                        }
                    }

                    tmpOutRe.Data[y, x, 0] = sumDataRe;
                    tmpOutIm.Data[y, x, 0] = sumDataIm;
                }
            }

            oRe = tmpOutRe.Clone();
            oIm = tmpOutIm.Clone();

            tmpOutRe.Dispose();
            tmpOutIm.Dispose();
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

        public void Convolution(ref Image<Gray, Single> input, ref Image<Gray, Single> output, ref Image<Gray, Single> mask)
        {
            Image<Gray, Single> tmpOut = new Image<Gray, float>(input.Width, input.Height);

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

                    tmpOut.Data[y, x, 0] = avgData;
                }
            }

            output = tmpOut.Clone();
            tmpOut.Dispose();
        }

        public void ScalingImage0To1(ref Image<Gray, Single> inout)
        {
            for (int x = 0; x < inout.Width; x++)
            {
                for (int y = 0; y < inout.Height; y++)
                {
                    if (!Double.IsNaN(inout.Data[y, x, 0]))
                    {
                        Gray gValue = new Gray();
                        gValue.Intensity = (inout.Data[y, x, 0] * 255.0);
                        inout[y, x] = gValue;
                    }
                    else
                    {
                        inout.Data[y, x, 0] = 0;
                    }
                }
            }
        }

        public void plotResponse(Image<Gray, Single> input, ref Image<Bgr, Single> output, int selectedLevel)
        {
            output = input.Convert<Bgr, Single>();
            output.ROI = new Rectangle((GaussianX.Width / 2), (GaussianX.Width / 2), input.Width - ((GaussianX.Width / 2) * 2), input.Height - ((GaussianX.Height / 2) * 2));
            output = output.Copy();
            Image<Bgr, Single> output2 = output.Clone();
            Bgr[] rectColor = new Bgr[6];
            rectColor[0] = new Bgr(Color.Red);
            rectColor[1] = new Bgr(Color.Green);
            rectColor[2] = new Bgr(Color.Blue);
            rectColor[3] = new Bgr(Color.MediumSpringGreen);
            rectColor[4] = new Bgr(Color.HotPink);
            rectColor[5] = new Bgr(Color.Black);
            for (int i = S1k.Length - 1, j = 0; i >= 0; i--, j++)
            {
                int multiplier = (int)Math.Pow(2, i);
                int offset = i * (S1k.Length - j) + 2;

                int x1 = S1kPMax[i][0].X * multiplier + offset;
                int y1 = S1kPMax[i][0].Y * multiplier + offset;

                int x2 = S2kPMax[i][0].X * multiplier + offset;
                int y2 = S2kPMax[i][0].Y * multiplier + offset;

                int x1W = S1kPMaxW[i][0].X * multiplier + offset;
                int y1W = S1kPMaxW[i][0].Y * multiplier + offset;

                int x2W = S2kPMaxW[i][0].X * multiplier + offset;
                int y2W = S2kPMaxW[i][0].Y * multiplier + offset;

                if (i == selectedLevel)
                {
                    float length = 10;
                    if (S1kMax[2][0] > 0.45)
                    {
                        coreW = new Point[1];
                        coreW[0] = new Point(x1W, y1W);
                        x1W = x1W - 6;
                        y1W = y1W - 6;
                        output.Draw(new Rectangle(x1W, y1W, 13, 13), rectColor[0], 1);

                        PointF str1 = new PointF(coreW[0].X, coreW[0].Y);
                        PointF end1 = new PointF();
                        end1.X = (float)(str1.X + length * Math.Cos(ThetaCore));
                        end1.Y = (float)(str1.Y + length * Math.Sin(ThetaCore));

                        LineSegment2DF CoreLine = new LineSegment2DF(str1, end1);
                        output.Draw(CoreLine, rectColor[4], 1);
                    }

                    if (S2kMax[2][0] > 0.50)
                    {
                        deltaW = new Point[1];
                        deltaW[0] = new Point(x2W, y2W);
                        output2.Draw(new Cross2DF(new PointF(x2W, y2W), 13, 13), rectColor[1], 1);

                        PointF str2 = new PointF(deltaW[0].X, deltaW[0].Y);
                        PointF end2 = new PointF();
                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta));

                        LineSegment2DF DeltaLine = new LineSegment2DF(str2, end2);
                        output2.Draw(DeltaLine, rectColor[3], 1);

                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta + (2 * Math.PI / 3)));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta + (2 * Math.PI / 3)));

                        DeltaLine = new LineSegment2DF(str2, end2);
                        output2.Draw(DeltaLine, rectColor[3], 1);

                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta - (2 * Math.PI / 3)));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta - (2 * Math.PI / 3)));

                        DeltaLine = new LineSegment2DF(str2, end2);
                        output2.Draw(DeltaLine, rectColor[3], 1);
                    }

                }
            }


            output = output.ConcateHorizontal(S1kDisp[0].Convert<Bgr, Single>());
            output = output.ConcateHorizontal(S1kDisp[1].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());
            output = output.ConcateHorizontal(S1kDisp[2].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());
            output = output.ConcateHorizontal(S1kDisp[3].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());

            output2 = output2.ConcateHorizontal(S2kDisp[0].Convert<Bgr, Single>());
            output2 = output2.ConcateHorizontal(S2kDisp[1].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());
            output2 = output2.ConcateHorizontal(S2kDisp[2].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());
            output2 = output2.ConcateHorizontal(S2kDisp[3].Resize(S1k[0].Width, S1k[0].Height, Inter.Nearest).Convert<Bgr, Single>());
            Image<Bgr, Single> output3 = output;
            output3 = output3.ConcateVertical(output2);

            output = output3.Clone();

            output2.Dispose();
            output3.Dispose();
        }

        public void plotResponseInOneImage(Image<Gray, Single> input, ref Image<Bgr, Single> output, int selectedLevel)
        {
            output = input.Convert<Bgr, Single>();
            output.ROI = new Rectangle((GaussianX.Width / 2), (GaussianX.Width / 2), input.Width - ((GaussianX.Width / 2) * 2), input.Height - ((GaussianX.Height / 2) * 2));
            output = output.Copy();
            Bgr[] rectColor = new Bgr[6];
            rectColor[0] = new Bgr(Color.Red);
            rectColor[1] = new Bgr(Color.Green);
            rectColor[2] = new Bgr(Color.Blue);
            rectColor[3] = new Bgr(Color.MediumSpringGreen);
            rectColor[4] = new Bgr(Color.HotPink);
            rectColor[5] = new Bgr(Color.Black);
            for (int i = S1k.Length - 1, j = 0; i >= 0; i--, j++)
            {
                int multiplier = (int)Math.Pow(2, i);
                int offset = i * (S1k.Length - j) + 2;
                int x1 = S1kPMax[i][0].X * multiplier + offset;
                int y1 = S1kPMax[i][0].Y * multiplier + offset;

                int x2 = S2kPMax[i][0].X * multiplier + offset;
                int y2 = S2kPMax[i][0].Y * multiplier + offset;

                int x1W = S1kPMaxW[i][0].X * multiplier + offset;
                int y1W = S1kPMaxW[i][0].Y * multiplier + offset;

                int x2W = S2kPMaxW[i][0].X * multiplier + offset;
                int y2W = S2kPMaxW[i][0].Y * multiplier + offset;

                if (i == selectedLevel)
                {
                    float length = 10;
                    if (S1kMax[2][0] > 0.45)
                    {
                        coreW = new Point[1];
                        coreW[0] = new Point(x1W, y1W);
                        x1W = x1W - 6;
                        y1W = y1W - 6;
                        output.Draw(new Rectangle(x1W, y1W, 13, 13), rectColor[0], 1);

                        PointF str1 = new PointF(coreW[0].X, coreW[0].Y);
                        PointF end1 = new PointF();
                        end1.X = (float)(str1.X + length * Math.Cos(ThetaCore));
                        end1.Y = (float)(str1.Y + length * Math.Sin(ThetaCore));

                        LineSegment2DF CoreLine = new LineSegment2DF(str1, end1);
                        output.Draw(CoreLine, rectColor[4], 1);
                    }
                    if (S2kMax[2][0] > 0.50)
                    {
                        deltaW = new Point[1];
                        deltaW[0] = new Point(x2W, y2W);
                        output.Draw(new Cross2DF(new PointF(x2W, y2W), 13, 13), rectColor[1], 1);

                        PointF str2 = new PointF(deltaW[0].X, deltaW[0].Y);
                        PointF end2 = new PointF();
                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta));

                        LineSegment2DF DeltaLine = new LineSegment2DF(str2, end2);
                        output.Draw(DeltaLine, rectColor[3], 1);

                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta + (2 * Math.PI / 3)));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta + (2 * Math.PI / 3)));

                        DeltaLine = new LineSegment2DF(str2, end2);
                        output.Draw(DeltaLine, rectColor[3], 1);

                        end2.X = (float)(str2.X + length * Math.Cos(ThetaDelta - (2 * Math.PI / 3)));
                        end2.Y = (float)(str2.Y + length * Math.Sin(ThetaDelta - (2 * Math.PI / 3)));

                        DeltaLine = new LineSegment2DF(str2, end2);
                        output.Draw(DeltaLine, rectColor[3], 1);
                    }
                }
            }
        }
    }
}
