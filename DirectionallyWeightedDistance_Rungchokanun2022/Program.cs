using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//For NaturalStringComparer
using System.Security;
using System.Runtime.InteropServices;
using System.IO;
//For PointF
using System.Drawing;
//For Image Utilities
using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.CV.CvEnum;
using Emgu.CV.UI;
//For Stopwatch
using System.Diagnostics;
//For Matched Minutiae Plot
using Graphs;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class Program
    {
        static void Main(string[] args)
        {
            //exampleSelectMinutiaeUnderPose();
            //exampleMinutiaeTripletsFormation();
            //exampleJYMatching();
            //exampleM3glMatching();
            //exampleLAS();
            //exampleComplexFilter();
            //exampleOrientationField();
        }

        public static void exampleMinutiaeTripletsFormation()
        {
            string FPImagePath = @"E:\Research\FVC2004_Db1FingerNet";//fingerprint image path
            string searchFPImagPattern = "*.bmp";//fingerprint image file extension
            string minutiaePath = FPImagePath + @"\FingerNet\Minutiae\"; // *.mnt files from FingerNet
            string ofPath = FPImagePath + @"\FingerNet\OF\"; // *.txt files which generated from FingerNetOFMatToText.m
            string fpSegmentPath = FPImagePath + @"\FingerNet\Segment\"; // *.png files from FingerNet

            string saveMTPath;//output path for minutiae-triplets text file and each file contains minutiae's id
            int Kneighbor;//number of neighbors to be selected to form minutiae-triplets
            int mode;//mode of minutiae-triplets formation
            bool arrangeCWandremoveDuplicate;//Set this flag to "true" to construct minutiae-triplaets for M3gl Algorithm, 
                                             // the minutiae in each minutiae-triplets are arranged in a clockwise direction 
                                             // and the duplicate minutiae-triplets are discarded.
                                             // Otherwise, Set this flag to "false".

            //uncomment to form minutiae-triplets using conventional 2-nearest neighbor selection method for JY matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\ConventionalNNMethod\";
            //Kneighbor = 2;
            //mode = 0;
            //arrangeCWandremoveDuplicate = false;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using conventional 2-nearest neighbor selection method for M3gl c=2 matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\ConventionalNNMethodM3glC2\";
            //Kneighbor = 2;
            //mode = 0;
            //arrangeCWandremoveDuplicate = true;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using directionally weighted distance 2-nearest neighbor selection method for JY matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\DirectionalWeightedNearestNeighbor_1stMethod\";
            //Kneighbor = 2;
            //mode = 1;
            //arrangeCWandremoveDuplicate = false;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using directionally weighted distance 2-nearest neighbor selection method for M3gl c=2 matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\DirectionalWeightedNearestNeighbor_1stMethodM3glC2\";
            //Kneighbor = 2;
            //mode = 1;
            //arrangeCWandremoveDuplicate = true;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using ridge flow directionally weighted distance 2-nearest neighbor selection method for JY matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\RidgeFlowDirectionalWeightedNearestNeighbor_2ndMethod2\";
            //Kneighbor = 2;
            //mode = 2;
            //arrangeCWandremoveDuplicate = false;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using ridge flow directionally weighted distance 2-nearest neighbor selection method for M3gl c=2 matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\RidgeFlowDirectionalWeightedNearestNeighbor_2ndMethodM3glC2\";
            //Kneighbor = 2;
            //mode = 2;
            //arrangeCWandremoveDuplicate = true;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using conventional 4-nearest neighbor selection method for M3gl c=4 matching algorithm (original M3gl c=4)
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\ConventionalNNMethodM3glC4\";
            //Kneighbor = 4;
            //mode = 0;
            //arrangeCWandremoveDuplicate = true;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using directionally weighted distance 4-nearest neighbor selection method for M3gl c=4 matching algorithm
            //saveMTPath = FPImagePath + @"\FingerNet\Minutiae\DirectionalWeightedNearestNeighbor_1stMethodM3glC4_34bw\";
            //Kneighbor = 4;
            //mode = 3;
            //arrangeCWandremoveDuplicate = true;
            //MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

            //uncomment to form minutiae-triplets using ridge flow directionally weighted distance 4-nearest neighbor selection method for M3gl c=4 matching algorithm
            saveMTPath = FPImagePath + @"\FingerNet\Minutiae\RidgeFlowDirectionalWeightedNearestNeighbor_2ndMethodM3glC4_34bw\";
            Kneighbor = 4;
            mode = 4;
            arrangeCWandremoveDuplicate = true;
            MinutiaeTripletsFormation(FPImagePath, searchFPImagPattern, minutiaePath, ofPath, fpSegmentPath, Kneighbor, mode, arrangeCWandremoveDuplicate, saveMTPath);

        }

        public static void exampleSelectMinutiaeUnderPose()
        {
            string FPImagePath = @"E:\Research\NIST14PNG_2700";//fingerprint image path
            string searchFPImagPattern = "*.png";//fingerprint image file extension
            string minutiaePath = FPImagePath + @"\FingerNet\Minutiae\"; // *.mnt files from FingerNet
            string posePath = FPImagePath + @"\Pose\"; // *.txt files which contains pose position (x,y) and direction (θ), 
            //where θ=[-90°,90°], 0° points in upward direction and increases in counter-clockwise direction.
            string saveMPath = FPImagePath + @"\FingerNet\MinutiaeUnderpose\"; // output path for minutiae .mnt file

            selectMinutiaeUnderPose(FPImagePath, searchFPImagPattern, minutiaePath, posePath, saveMPath);
        }

        public static void exampleJYMatching()
        {
            string DbConsole = "FVC2004_Db1";
            string FPImagePath = @"E:\Research\FVC2004_Db1FingerNet";//fingerprint image path
            string fpImgExtension = ".bmp";//fingerprint image file extension
            string minutiaePath = FPImagePath + @"\FingerNet\Minutiae\"; // *.mnt files from FingerNet
            string minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\ConventionalNNMethod\"; // minutiae-triplets *.txt files from MinutiaeTripletsFormation()
            string saveResultPath = FPImagePath + @"\FingerNet\ResultToGitHub\"; //output path to store matching results
            string saveResultFilename = "JYFNMRNoRCMT_onConventionalNNMethod"; //matching results filename which is saved as .txt
            int numOfImagePerFinger = 8;
            //JY
            Console.WriteLine("Running JY on " + DbConsole);
            matchJYFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, true);
            saveResultFilename = "JYFMRNoRCMT_onConventionalNNMethod";
            matchJYFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

            minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\DirectionalWeightedNearestNeighbor_1stMethod\";
            saveResultFilename = "JYFNMRNoRCMT_on1stMethod";
            matchJYFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, false);
            saveResultFilename = "JYFMRNoRCMT_on1stMethod";
            matchJYFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

            minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\RidgeFlowDirectionalWeightedNearestNeighbor_2ndMethod2\";
            saveResultFilename = "JYFNMRNoRCMT_on2ndMethod";
            matchJYFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, false);
            saveResultFilename = "JYFMRNoRCMT_on2ndMethod";
            matchJYFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

        }

        public static void exampleM3glMatching()
        {
            string DbConsole = "FVC2004_Db1";
            string FPImagePath = @"E:\Research\FVC2004_Db1FingerNet";//fingerprint image path
            string fpImgExtension = ".bmp";//fingerprint image file extension
            string minutiaePath = FPImagePath + @"\FingerNet\Minutiae\"; // *.mnt files from FingerNet
            string minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\ConventionalNNMethodM3glC2\"; // minutiae-triplets *.txt files from MinutiaeTripletsFormation()
            string saveResultPath = FPImagePath + @"\FingerNet\ResultToGitHub\"; //output path to store matching results
            string saveResultFilename = "M3glFNMR_onConventionalNNMethodM3glC2"; //matching results filename which is saved as .txt
            int numOfImagePerFinger = 8;
            //M3gl
            Console.WriteLine("Running M3gl on " + DbConsole);
            matchM3glFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, false);
            saveResultFilename = "M3glFMR_onConventionalNNMethodM3glC2";
            matchM3glFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

            minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\DirectionalWeightedNearestNeighbor_1stMethodM3glC2\";
            saveResultFilename = "M3glFNMR_on1stMethodC2";
            matchM3glFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, false);
            saveResultFilename = "M3glFMR_on1stMethodC2";
            matchM3glFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

            minutiaeTripletsPath = FPImagePath + @"\FingerNet\Minutiae\RidgeFlowDirectionalWeightedNearestNeighbor_2ndMethodM3glC2\";
            saveResultFilename = "M3glFNMR_on2ndMethodC2";
            matchM3glFNMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, numOfImagePerFinger, false);
            saveResultFilename = "M3glFMR_on2ndMethodC2";
            matchM3glFMR(minutiaePath, minutiaeTripletsPath, FPImagePath, fpImgExtension, saveResultPath, saveResultFilename, false);

        }

        public static void exampleLAS()
        {
            LAS_Liu2006 las = new LAS_Liu2006();
            las.runLAS(@"E:\Research\FVC2000_Db2\1_1.tif", @"‪E:\Research\FVC2000_Db2\LAS");
        }

        public static void exampleComplexFilter()
        {
            Image<Gray, Single> input = new Image<Gray, Single>(@"E:\Research\FVC2000_Db2\41_4.tif");
            Image<Bgr, Single> output = new Image<Bgr, float>(1, 1);
            ComplexFilter_Nilson2003 cpFilter = new ComplexFilter_Nilson2003();
            cpFilter.Initialize();
            cpFilter.calFilterResponse(input, 20);
            cpFilter.plotResponse(input, ref output, 1);
            ImageViewer.Show(output, "Filter Response");

            cpFilter.free();
            output.Dispose();
            input.Dispose();
        }

        public static void exampleOrientationField()
        {
            Image<Gray, Single> input = new Image<Gray, Single>(@"E:\Research\FVC2000_Db2\41_4.tif");
            Image<Gray, Single> Theta = new Image<Gray, float>(1, 1);
            Image<Gray, Single> Coh = new Image<Gray, float>(1, 1);
            OF_Bazen2002 ofBazen = new OF_Bazen2002();

            ofBazen.calOF(input, ref Theta, ref Coh, 20);

            Image<Bgr, Byte> ThetaJet = Theta.Convert<Bgr, Byte>();
            CvInvoke.ApplyColorMap(ThetaJet, ThetaJet, ColorMapType.Hsv);
            ImageViewer.Show(ThetaJet);

            int blockSize = 16;
            Image<Gray, Single> ofDownImg = new Image<Gray, float>(1, 1);
            Image<Gray, Single> cohDownImg = new Image<Gray, float>(1, 1); 
            Image<Gray, Single> OutputQuiver = new Image<Gray, Single>(1, 1);
            ofBazen.averageOF(ref Theta, ref Coh, ref ofDownImg, ref cohDownImg, blockSize);
            Quiver(ref ofDownImg, ref cohDownImg, ref OutputQuiver, 15, false);
            ImageViewer.Show(OutputQuiver.Convert<Bgr, Single>());

            OutputQuiver.Dispose();
            cohDownImg.Dispose();
            ofDownImg.Dispose();
            
            ThetaJet.Dispose();
            
            ofBazen.free();
            
            Coh.Dispose();
            Theta.Dispose();
            input.Dispose();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="fpImgPath"></param>
        /// <param name="searchFPImagPattern"></param>
        /// <param name="minutiaePath"></param>
        /// <param name="ofPath"></param>
        /// <param name="fpSegmentPath"></param>
        /// <param name="Kneighbor">number of K nearest neighbor to be selected and used to construct minutiae-triplets</param>
        /// <param name="mode"> "0" = conventional k-nearest neighbor
        ///                     "1" = directionally weighted distance for 2-nearest neighbor selection, 
        ///                           set arrangeCWandremoveDuplicate to false for JY algorithm
        ///                           set arrangeCWandremoveDuplicate to true for M3gl algorithm
        ///                     "2" = ridge flow directionally weighted distance for 2-nearest neighbor selection,
        ///                           set arrangeCWandremoveDuplicate to false for JY algorithm
        ///                           set arrangeCWandremoveDuplicate to true for M3gl algorithm
        ///                     "3" = directionally weighted distance for 4-nearest neighbor selection,
        ///                           alway set arrangeCWandremoveDuplicate to true for M3gl algorithm
        ///                     "4" = ridge flow directionally weighted distance for 4-nearest neighbor selection,
        ///                           alway set arrangeCWandremoveDuplicate to true for M3gl algorithm
        /// </param>
        /// <param name="arrangeCWandremoveDuplicate">
        /// Set this flag to "true" to construct minutiae-triplaets for M3gl Algorithm, 
        /// the minutiae in each minutiae-triplets are arranged in a clockwise direction 
        /// and the duplicate minutiae-triplets are discarded.
        /// Otherwise, Set this flag to "false".</param>
        /// <param name="saveMPath"></param>
        public static void MinutiaeTripletsFormation(string fpImgPath, string searchFPImagPattern, string minutiaePath, string ofPath, string fpSegmentPath,
                                        int Kneighbor, int mode, bool arrangeCWandremoveDuplicate, string saveMTPath)
        {
            var vFPFileNames = new DirectoryInfo(fpImgPath).GetFileSystemInfos(searchFPImagPattern).OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);

            string[] fpFileNames = vFPFileNames.ToArray();

            double D = 84.0;
            double Theta = Math.PI / 2.0;
            double L = 84.0;

            int ofBlockSize = 8;

            for (int f = 0; f < fpFileNames.Length; f++)
            {
                string fname = Path.GetFileNameWithoutExtension(fpFileNames[f]);

                Image<Bgr, Single> inputBgr = new Image<Bgr, Single>(fpFileNames[f]);

                Size fpImgSize = new Size();
                KMinutia[] M = readFingerNetMinutiaeAndPlot(minutiaePath + fname + ".mnt", ref fpImgSize, ref inputBgr, false);

                List<mTriplet> listMTriplets = new List<mTriplet>();
                switch (mode) {
                    case 0:
                        if (arrangeCWandremoveDuplicate)
                        {
                            listMTriplets = conventionalNearestNeighborM3glArrangeClockwise(M, Kneighbor);
                        }
                        else
                        {
                            listMTriplets = conventionalNearestNeighbor(M, Kneighbor);
                        }
                        break;
                    case 1:
                        listMTriplets = directionallyWeightedDistanceNearestNeighbor_Method1(M, D, Theta, Kneighbor, arrangeCWandremoveDuplicate);
                        break;
                    case 2:
                        {
                            Image<Gray, Single> inputSeg = new Image<Gray, float>(fpSegmentPath + fname + "_seg.png");
                            Image<Gray, Single> inputSegBlock = new Image<Gray, float>(1, 1);

                            pixelToBlock(ref inputSeg, ref inputSegBlock, ofBlockSize);

                            Image<Gray, Single> ofImgBlock = new Image<Gray, Single>(1, 1);

                            getOFFromFile(ofPath + fname + ".txt", inputBgr.Size, ref ofImgBlock, ofBlockSize, inputSegBlock);

                            #region plot OF Quiver
                            //Image<Gray, Single> ofImg = new Image<Gray, Single>(1, 1);
                            //BlocksToPixels(ref ofImgBlock, ref ofImg, ofBlockSize);
                            //Image<Gray, Single> cohImgBlock = new Image<Gray, Single>(ofImgBlock.Width, ofImgBlock.Height, new Gray(1));
                            //Image<Gray, Single> OutputQuiverTheta = new Image<Gray, Single>(1, 1);

                            //Quiver(ref ofImgBlock, ref cohImgBlock, ref OutputQuiverTheta, ofBlockSize);

                            //double alpha = 0.5;
                            //double beta = (1.0 - alpha);
                            //double gamma = 0.0;
                            //inputBgr = inputBgr.AddWeighted(OutputQuiverTheta.Convert<Bgr, Single>(), alpha, beta, gamma);
                            //ImageViewer.Show(inputBgr.ConcateHorizontal(OutputQuiverTheta.Convert<Bgr, Single>()));
                            
                            //OutputQuiverTheta.Dispose();
                            //cohImgBlock.Dispose();
                            //ofImg.Dispose();
                            #endregion

                            listMTriplets = ridgeFlowDirectionallyWeightedDistanceNearestNeighbor_Method2(M, ofImgBlock, ofBlockSize, D, Theta, L, Kneighbor, arrangeCWandremoveDuplicate);

                            inputSeg.Dispose();
                            inputSegBlock.Dispose();
                            ofImgBlock.Dispose();
                        }
                        break;
                    case 3:
                        listMTriplets = directionallyWeightedDistanceNearestNeighbor_Method1_34bw(M, D, Theta, Kneighbor, arrangeCWandremoveDuplicate);
                        break;
                    case 4:
                        {
                            Image<Gray, Single> inputSeg = new Image<Gray, float>(fpSegmentPath + fname + "_seg.png");
                            Image<Gray, Single> inputSegBlock = new Image<Gray, float>(1, 1);

                            pixelToBlock(ref inputSeg, ref inputSegBlock, ofBlockSize);

                            Image<Gray, Single> ofImgBlock = new Image<Gray, Single>(1, 1);

                            getOFFromFile(ofPath + fname + ".txt", inputBgr.Size, ref ofImgBlock, ofBlockSize, inputSegBlock);

                            #region plot OF Quiver
                            //Image<Gray, Single> ofImg = new Image<Gray, Single>(1, 1);
                            //BlocksToPixels(ref ofImgBlock, ref ofImg, ofBlockSize);
                            //Image<Gray, Single> cohImgBlock = new Image<Gray, Single>(ofImgBlock.Width, ofImgBlock.Height, new Gray(1));
                            //Image<Gray, Single> OutputQuiverTheta = new Image<Gray, Single>(1, 1);

                            //Quiver(ref ofImgBlock, ref cohImgBlock, ref OutputQuiverTheta, ofBlockSize);

                            //double alpha = 0.5;
                            //double beta = (1.0 - alpha);
                            //double gamma = 0.0;
                            //inputBgr = inputBgr.AddWeighted(OutputQuiverTheta.Convert<Bgr, Single>(), alpha, beta, gamma);
                            //ImageViewer.Show(inputBgr.ConcateHorizontal(OutputQuiverTheta.Convert<Bgr, Single>()));

                            //OutputQuiverTheta.Dispose();
                            //cohImgBlock.Dispose();
                            //ofImg.Dispose();
                            #endregion

                            listMTriplets = ridgeFlowDirectionallyWeightedDistanceNearestNeighbor_Method2_34bw(M, ofImgBlock, ofBlockSize, D, Theta, L, Kneighbor, arrangeCWandremoveDuplicate);
                            
                            inputSeg.Dispose();
                            inputSegBlock.Dispose();
                            ofImgBlock.Dispose();
                        }
                        break;
                    default: break;
                }

                saveMTripletsToText(listMTriplets, saveMTPath + fname + ".txt");

                listMTriplets.Clear();
                listMTriplets.TrimExcess();
                listMTriplets = null;
                inputBgr.Dispose();
            }
        }

        public static void selectMinutiaeUnderPose(string fpImgPath, string searchFPImagPattern, string minutiaePath, string posePath, string saveMPath)
        {
            var vFPFileNames = new DirectoryInfo(fpImgPath).GetFileSystemInfos(searchFPImagPattern).Select(x => x.FullName);
            string[] fpFileNames = vFPFileNames.ToArray();

            for (int f = 0; f < fpFileNames.Length; f++)
            {
                string fname = Path.GetFileNameWithoutExtension(fpFileNames[f]);

                Image<Bgr, Single> inputBgr = new Image<Bgr, Single>(fpFileNames[f]);
                Size size = new Size();

                KMinutia[] M = readFingerNetMinutiaeAndPlot(minutiaePath + fname + ".mnt", ref size, ref inputBgr, false);

                Pose pose = getPoseFromFile(posePath + fname + ".txt", ref inputBgr, 14);

                determineMinutiaeLowerOrUpperFromPoseBase(M, pose, true);
                var mUnderPose = M.Where(m => m.valid).ToList();
                mUnderPose.Sort((a, b) => a.M_ID.CompareTo(b.M_ID));

                #region draw pose and its base on FPImage
                //drawPoseOnImage(ref inputBgr, pose);
                //PointF endP1 = new PointF();
                //PointF endP2 = new PointF();
                //double poseLength = 500;
                //endP1.X = (float)(pose.p.X + poseLength * Math.Cos(pose.angle + (Math.PI / 2.0)));
                //endP1.Y = (float)(pose.p.Y + poseLength * Math.Sin(pose.angle + (Math.PI / 2.0)));
                //endP2.X = (float)(pose.p.X + poseLength * Math.Cos(pose.angle - (Math.PI / 2.0)));
                //endP2.Y = (float)(pose.p.Y + poseLength * Math.Sin(pose.angle - (Math.PI / 2.0)));

                //LineSegment2DF poseBase = new LineSegment2DF(endP1, endP2);
                //inputBgr.Draw(poseBase, new Bgr(Color.Red), 2);
                //int i = 0;
                //foreach (KMinutia m in mUnderPose)
                //{
                //    PointF str = new PointF(m.p.X, m.p.Y);
                //    PointF end = new PointF();
                //    double length = 15;
                //    end.X = (float)(str.X + length * Math.Cos(m.direction));
                //    end.Y = (float)(str.Y + length * Math.Sin(m.direction));
                //    LineSegment2DF minutiaLine = new LineSegment2DF(str, end);
                //    inputBgr.Draw(new CircleF(m.p, 5), new Bgr(Color.Red), 1);
                //    inputBgr.Draw(minutiaLine, new Bgr(Color.Red), 1);
                //    //inputBgr.Draw(m.M_ID.ToString(), new Point((int)str.X - 5, (int)str.Y + 13), FontFace.HersheyPlain, 0.5, new Bgr(Color.Red));
                //    inputBgr.Draw(i.ToString(), new Point((int)str.X - 5, (int)str.Y + 13), FontFace.HersheyPlain, 0.5, new Bgr(Color.Red));
                //    i++;
                //}
                //ImageViewer.Show(inputBgr);
                #endregion

                string savePath = saveMPath + fname + ".mnt";

                using (var writer = new StreamWriter(savePath))
                {
                    writer.WriteLine("{0}", fname);
                    writer.WriteLine("{0} {1} {2}", mUnderPose.Count.ToString(), inputBgr.Height.ToString(), inputBgr.Width.ToString());

                    for (int m = 0; m < mUnderPose.Count; m++)
                    {
                        double angle = mUnderPose[m].direction;
                        writer.WriteLine("{0} {1} {2}", mUnderPose[m].p.X, mUnderPose[m].p.Y, angle);
                    }
                }

                mUnderPose.Clear();
                mUnderPose.TrimExcess();
                mUnderPose = null;
                inputBgr.Dispose();
            }
        }

        public static KMinutia[] readFingerNetMinutiaeAndPlot(string mPath, ref Size size, ref Image<Bgr, Single> inputBgr, bool plotMinuntiae = false)
        {
            string[] strMinutiaLines = System.IO.File.ReadAllLines(mPath);

            char[] splitChars = { '\t', ' ' };

            string fname = strMinutiaLines[0];
            string[] strResNNumOfMinutiae = strMinutiaLines[1].Split(splitChars);

            int numOfMinutiae = 0;

            int imgWidth = Convert.ToInt32(strResNNumOfMinutiae[2]);
            int imgHeight = Convert.ToInt32(strResNNumOfMinutiae[1]);
            size.Width = imgWidth;
            size.Height = imgHeight;

            if (Int32.TryParse(strResNNumOfMinutiae[0], out numOfMinutiae) && numOfMinutiae != 0)
            {
                KMinutia[] M = new KMinutia[numOfMinutiae];

                for (int l = 2, cM = 0; l < strMinutiaLines.Length; l++, cM++)
                {
                    string[] strMinutia = strMinutiaLines[l].Split(splitChars);

                    int x = Convert.ToInt32(strMinutia[0]);
                    int y = Convert.ToInt32(strMinutia[1]);
                    double direction = Convert.ToDouble(strMinutia[2]);

                    M[cM] = new KMinutia();
                    M[cM].M_ID = cM;
                    M[cM].p = new Point(x, y);
                    M[cM].direction = direction;

                    #region plot original minutiae
                    if (plotMinuntiae)
                    {
                        PointF str = new PointF(M[cM].p.X, M[cM].p.Y);
                        PointF end = new PointF();
                        double length = 15;
                        end.X = (float)(str.X + length * Math.Cos(direction));
                        end.Y = (float)(str.Y + length * Math.Sin(direction));
                        LineSegment2DF minutiaLine = new LineSegment2DF(str, end);
                        inputBgr.Draw(new CircleF(M[cM].p, 5), new Bgr(0, 0, 255), 1);
                        inputBgr.Draw(minutiaLine, new Bgr(0, 0, 255), 1);
                        inputBgr.Draw(M[cM].M_ID.ToString(), new Point((int)str.X - 5, (int)str.Y + 13), FontFace.HersheyPlain, 0.5, new Bgr(0, 0, 255));

                        ////Arrow minutiae plot
                        //length = 20;
                        //end.X = (float)(str.X + length * Math.Cos(direction));
                        //end.Y = (float)(str.Y + length * Math.Sin(direction));
                        //minutiaLine = new LineSegment2DF(str, end);
                        //double lengthArrow = 7;
                        //PointF endArrow1 = new PointF();
                        //endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(direction + (Math.PI / 7)));
                        //endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(direction + (Math.PI / 7)));
                        //PointF endArrow2 = new PointF();
                        //endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(direction - (Math.PI / 7)));
                        //endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(direction - (Math.PI / 7)));
                        //LineSegment2DF arrowLine1 = new LineSegment2DF(end, endArrow1);
                        //LineSegment2DF arrowLine2 = new LineSegment2DF(end, endArrow2);
                        //inputBgr.Draw(new CircleF(M[cM].p, 5), new Bgr(255, 0, 0), 2);
                        //inputBgr.Draw(minutiaLine, new Bgr(255, 0, 0), 2);
                        //inputBgr.Draw(arrowLine1, new Bgr(255, 0, 0), 2);
                        //inputBgr.Draw(arrowLine2, new Bgr(255, 0, 0), 2);
                    }
                    #endregion
                }

                return M;
            }

            return new KMinutia[numOfMinutiae];
        }

        public static KMinutia[] readFingerNetMinutiae(string mPath)
        {
            string[] strMinutiaLines = System.IO.File.ReadAllLines(mPath);

            char[] splitChars = { '\t', ' ' };

            string fname = strMinutiaLines[0];
            string[] strResNNumOfMinutiae = strMinutiaLines[1].Split(splitChars);

            int numOfMinutiae = 0;

            int imgWidth = Convert.ToInt32(strResNNumOfMinutiae[2]);
            int imgHeight = Convert.ToInt32(strResNNumOfMinutiae[1]);

            if (Int32.TryParse(strResNNumOfMinutiae[0], out numOfMinutiae) && numOfMinutiae != 0)
            {
                KMinutia[] M = new KMinutia[numOfMinutiae];

                for (int l = 2, cM = 0; l < strMinutiaLines.Length; l++, cM++)
                {
                    string[] strMinutia = strMinutiaLines[l].Split(splitChars);

                    int x = Convert.ToInt32(strMinutia[0]);
                    int y = Convert.ToInt32(strMinutia[1]);
                    double direction = Convert.ToDouble(strMinutia[2]);

                    M[cM] = new KMinutia();
                    M[cM].M_ID = cM;
                    M[cM].p = new Point(x, y);
                    M[cM].direction = direction;
                }

                return M;
            }

            return new KMinutia[numOfMinutiae];
        }

        public static List<mTriplet> readMinutiaeTripletsFromTextFile(ref KMinutia[] M, string mTripletPath)
        {
            List<int> inputM = M.Select(m => m.M_ID).ToList();
            List<KMinutia> validM = M.ToList();

            List<mTriplet> T = new List<mTriplet>();

            string[] pID_mTripletLines = System.IO.File.ReadAllLines(mTripletPath);

            char[] splitChars = { '\t', ' ' };

            for (int l = 0; l < pID_mTripletLines.Length; l++)
            {
                string[] pID_mTriplet = pID_mTripletLines[l].Split(splitChars);

                int p0 = Convert.ToInt32(pID_mTriplet[0]);
                int p1 = Convert.ToInt32(pID_mTriplet[1]);
                int p2 = Convert.ToInt32(pID_mTriplet[2]);

                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(M, p => p.M_ID == p0);
                int idxP1 = Array.FindIndex(M, p => p.M_ID == p1);
                int idxP2 = Array.FindIndex(M, p => p.M_ID == p2);
                t.pi[0] = M[idxP0];
                t.pi[1] = M[idxP1];
                t.pi[2] = M[idxP2];

                T.Add(t);

                inputM.Remove(p0);
                inputM.Remove(p1);
                inputM.Remove(p2);
            }

            for (int i = 0; i < inputM.Count; i++)
            {
                validM.RemoveAll(m => m.M_ID == inputM[i]);
            }
            M = validM.ToArray();

            inputM.Clear();
            inputM.TrimExcess();
            inputM = null;

            return T;
        }

        public static Pose getPoseFromFile(string path, ref Image<Bgr, Single> inputBgr, int NISTDb)
        {
            Pose p = new Pose();
            string[] strPoseLines = System.IO.File.ReadAllLines(path);
            char[] splitChars = { '\t', ' ' };
            string[] stringPose = strPoseLines[0].Split(splitChars);
            p.p.X = Convert.ToInt32(stringPose[0]);
            p.p.Y = Convert.ToInt32(stringPose[1]);

            if (NISTDb == 14)
            {
                p.angle = (-0.5 * Math.PI) - (Convert.ToDouble(stringPose[2]) * Math.PI / 180.0);
            }
            else
            {
                p.angle = Convert.ToDouble(stringPose[2]);
            }

            return p;
        }

        public static void drawPoseOnImage(ref Image<Bgr, Single> img, Pose pose)
        {
            PointF strPose = new PointF(pose.p.X, pose.p.Y);
            PointF endPose = new PointF();
            PointF endArrowPose1 = new PointF();
            PointF endArrowPose2 = new PointF();
            double lengthPose = 200;
            double lengthArrowPose = 30;
            endPose.X = (float)(strPose.X + lengthPose * Math.Cos(pose.angle));
            endPose.Y = (float)(strPose.Y + lengthPose * Math.Sin(pose.angle));
            endArrowPose1.X = (float)(endPose.X - lengthArrowPose * Math.Cos(pose.angle + Math.PI / 7));
            endArrowPose1.Y = (float)(endPose.Y - lengthArrowPose * Math.Sin(pose.angle + Math.PI / 7));
            endArrowPose2.X = (float)(endPose.X - lengthArrowPose * Math.Cos(pose.angle - Math.PI / 7));
            endArrowPose2.Y = (float)(endPose.Y - lengthArrowPose * Math.Sin(pose.angle - Math.PI / 7));

            LineSegment2DF poseLine1 = new LineSegment2DF(strPose, endPose);
            LineSegment2DF arrowLine11 = new LineSegment2DF(endPose, endArrowPose1);
            LineSegment2DF arrowLine12 = new LineSegment2DF(endPose, endArrowPose2);

            img.Draw(new CircleF(pose.p, 5), new Bgr(0, 255, 255), 0);
            img.Draw(poseLine1, new Bgr(0, 255, 255), 3);
            img.Draw(arrowLine11, new Bgr(0, 255, 255), 3);
            img.Draw(arrowLine12, new Bgr(0, 255, 255), 3);
        }

        /// <summary>
        /// Determine that each minutia is lower or upper of pose's base.
        /// </summary>
        /// <param name="M"></param>
        /// <param name="pose"></param>
        /// <param name="lowerUpper">true = minutiae lying lower the pose base are valid, false = minutiae lying upper the pose base are valid</param>
        public static void determineMinutiaeLowerOrUpperFromPoseBase(KMinutia[] M, Pose pose, bool lowerUpper)
        {
            PointF endP1 = new PointF();
            PointF endP2 = new PointF();
            double poseLength = 500;
            endP1.X = (float)(pose.p.X + poseLength * Math.Cos(pose.angle + (Math.PI / 2.0)));
            endP1.Y = (float)(pose.p.Y + poseLength * Math.Sin(pose.angle + (Math.PI / 2.0)));
            endP2.X = (float)(pose.p.X + poseLength * Math.Cos(pose.angle - (Math.PI / 2.0)));
            endP2.Y = (float)(pose.p.Y + poseLength * Math.Sin(pose.angle - (Math.PI / 2.0)));

            for (int m = 0; m < M.Length; m++)
            {
                M[m].valid = determineLeftOrRightOfLine(endP1, endP2, M[m].p, lowerUpper);
            }
        }

        public static bool determineLeftOrRightOfLine(PointF lineP1, PointF lineP2, PointF queryP, bool lowerUpper)
        {
            double result = ((lineP2.X - lineP1.X) * (queryP.Y - lineP1.Y) - (lineP2.Y - lineP1.Y) * (queryP.X - lineP1.X));

            if (lowerUpper)
            {
                return result < 0 ? true : false;
            }
            else
            {
                return result > 0 ? true : false;
            }
        }

        public static List<mTriplet> conventionalNearestNeighbor(KMinutia[] P, int n = 2)
        {
            calDistanceOfPiToAllP(P);

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> mT = new List<mTriplet>();
            int mPerTriplet = 3;

            if (P.Length >= mPerTriplet)
            {
                foreach (KMinutia pi in P)
                {
                    var pID = pi.distanceFromMList.Keys.Take(n).ToList();
                    pID.Insert(0, pi.M_ID);
                    pID_mTriplet.Add(pID.ToArray());
                }

                for (int i = 0; i < pID_mTriplet.Count; i++)
                {
                    mTriplet t = new mTriplet();
                    int idxP0 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][0]);
                    int idxP1 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][1]);
                    int idxP2 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][2]);

                    t.pi[0] = P[idxP0];
                    t.pi[1] = P[idxP1];
                    t.pi[2] = P[idxP2];
                    mT.Add(t);
                }
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            return mT;
        }

        public static List<mTriplet> conventionalNearestNeighborM3glArrangeClockwise(KMinutia[] P, int n = 4)
        {
            calDistanceOfPiToAllP(P);

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> T = new List<mTriplet>();
            int mPerTriplet = 3;

            foreach (KMinutia pi in P)
            {
                var pID = pi.distanceFromMList.Keys.Take(n).ToList();
                pID.Insert(0, pi.M_ID);
                if (pID.Count >= mPerTriplet)
                {
                    int i = 0;
                    for (int j = 1; j < pID.Count; j++)
                    {
                        for (int k = j + 1; k < pID.Count; k++)
                        {
                            int[] inmID = new int[3];
                            inmID[0] = pID[i];
                            inmID[1] = pID[j];
                            inmID[2] = pID[k];

                            int[] outmID = arrangeClockwise(P, inmID);

                            pID_mTriplet.Add(outmID);
                        }
                    }
                }
            }

            pID_mTriplet = pID_mTriplet.Distinct(new IntArrayComparerIndexSensitive()).ToList();

            for (int i = 0; i < pID_mTriplet.Count; i++)
            {
                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][0]);
                int idxP1 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][1]);
                int idxP2 = Array.FindIndex(P, p => p.M_ID == pID_mTriplet[i][2]);

                t.pi[0] = P[idxP0];
                t.pi[1] = P[idxP1];
                t.pi[2] = P[idxP2];
                T.Add(t);
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            return T;
        }

        public static List<mTriplet> directionallyWeightedDistanceNearestNeighbor_Method1(KMinutia[] M, double D, double Theta, int n = 2, bool arrangeCWandremoveDuplicate = false)
        {
            for (int cM = 0; cM < M.Length; cM++)
            {
                M[cM].distanceFromMList = new Dictionary<int, double>();
                M[cM].radialAngleFromMList = new Dictionary<int, double>();

                for (int iM = 0; iM < M.Length; iM++)
                {
                    if (cM != iM)
                    {
                        PointF p1 = M[cM].p;
                        double directionP1 = M[cM].direction;

                        PointF p2 = M[iM].p;

                        double distance = distanceOf2Points(p1, p2);
                        double polarAngle = angleOf2Points(p2, p1);
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, directionP1);
                        M[cM].distanceFromMList.Add(M[iM].M_ID, distance);
                        M[cM].radialAngleFromMList.Add(M[iM].M_ID, radialAngle);
                    }
                }
                M[cM].distanceFromMList = M[cM].distanceFromMList.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
            }

            List<int> inputM = M.Select(m => m.M_ID).ToList();

            int mPerTriplet = 3;

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> T = new List<mTriplet>();

            if (inputM.Count >= mPerTriplet)
            {
                for (int i = 0; i < inputM.Count; i++)
                {
                    int midx = Array.FindIndex(M, m => m.M_ID == inputM[i]);

                    Dictionary<int, double> distance = new Dictionary<int, double>();
                    distance = M[inputM[i]].distanceFromMList.ToDictionary(k => k.Key, v => v.Value);
                    distance = distance.OrderBy(u => u.Value).ToDictionary(k => k.Key, v => v.Value);
                    var key = distance.Keys.ToArray();

                    var id = distance.Keys.ToArray();
                    var distanceList = distance.Values.ToArray();

                    double[] wRadius = new double[id.Length];
                    double[] wTheta = new double[id.Length];
                    double[] wScore = new double[id.Length];

                    for (int k = 0; k < id.Length; k++)
                    {
                        wRadius[k] = (distanceList[k] / D);
                        wTheta[k] = M[midx].radialAngleFromMList[id[k]] / Theta;
                        wScore[k] = Math.Sqrt(Math.Pow(wRadius[k], 2) + Math.Pow(wTheta[k], 2));
                    }
                    Array.Sort(wScore, id);

                    var d = id.Take(n).ToList();

                    d.Insert(0, inputM[i]);

                    if (arrangeCWandremoveDuplicate)
                    {
                        #region arrangeClockwise M3gl
                        //M3gl
                        //List<int[]> tmp_mTriplets = new List<int[]>();
                        //tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(pID.ToArray(), mPerTriplet)); //same result as below for-loop
                        int ii = 0;
                        for (int j = 1; j < d.Count; j++)
                        {
                            for (int k = j + 1; k < d.Count; k++)
                            {
                                int[] inmID = new int[3];
                                inmID[0] = d[ii];
                                inmID[1] = d[j];
                                inmID[2] = d[k];

                                int[] outmID = arrangeClockwise(M, inmID);

                                pID_mTriplet.Add(outmID);
                            }
                        }
                        #endregion
                    }
                    else
                    {
                        #region order the first neighboring minutia and the second neighboring minutia based on Jea and Govindaraju (2005)
                        List<int[]> tmp_mTriplets = new List<int[]>();
                        tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(d.ToArray(), mPerTriplet));
                        for (int mt = 0; mt < tmp_mTriplets.Count; mt++)
                        {
                            int[] pID = new int[3];
                            int idxP0 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][0]);
                            int idxP1 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][1]);
                            int idxP2 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][2]);

                            //AR
                            double signAcuteAngle = crossProductRightHandRule(M[idxP0].p, M[idxP1].p, M[idxP2].p);
                            if (signAcuteAngle >= 0)
                            {
                                pID[0] = M[idxP0].M_ID;
                                pID[1] = M[idxP1].M_ID;
                                pID[2] = M[idxP2].M_ID;
                            }
                            else
                            {
                                pID[0] = M[idxP0].M_ID;
                                pID[1] = M[idxP2].M_ID;
                                pID[2] = M[idxP1].M_ID;
                            }

                            pID_mTriplet.Add(pID);
                        }
                        tmp_mTriplets.Clear();
                        tmp_mTriplets.TrimExcess();
                        tmp_mTriplets = null;
                        #endregion
                    }

                    d.Clear();
                    d.TrimExcess();
                    d = null;

                    distance.Clear();
                    distance = null;
                }
            }

            if (arrangeCWandremoveDuplicate)
            {
                pID_mTriplet = pID_mTriplet.Distinct(new IntArrayComparerIndexSensitive()).ToList();
            }

            for (int i = 0; i < pID_mTriplet.Count; i++)
            {
                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][0]);
                int idxP1 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][1]);
                int idxP2 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][2]);

                t.pi[0] = M[idxP0];
                t.pi[1] = M[idxP1];
                t.pi[2] = M[idxP2];

                //var parentList = pID_mTriplet.Where(m => m[1] == t.pi[0].M_ID || m[2] == t.pi[0].M_ID).Select(mid => mid[0]).ToList();

                //t.parent = parentList.Distinct().ToList();

                T.Add(t);
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            inputM.Clear();
            inputM.TrimExcess();
            inputM = null;

            return T;
        }

        public static List<mTriplet> ridgeFlowDirectionallyWeightedDistanceNearestNeighbor_Method2(KMinutia[] M, Image<Gray, Single> ofImgBlock, int ofBlockSize, double D, double Theta, double L, int n = 2, bool arrangeCWandremoveDuplicate = false)
        {
            for (int cM = 0; cM < M.Length; cM++)
            {
                M[cM].distanceFromMList = new Dictionary<int, double>();
                M[cM].radialAngleFromMList = new Dictionary<int, double>();

                for (int iM = 0; iM < M.Length; iM++)
                {
                    if (cM != iM)
                    {
                        PointF p1 = M[cM].p;
                        double directionP1 = M[cM].direction;

                        PointF p2 = M[iM].p;

                        double distance = distanceOf2Points(p1, p2);
                        double polarAngle = angleOf2Points(p2, p1);
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, directionP1);
                        M[cM].distanceFromMList.Add(M[iM].M_ID, distance);
                        M[cM].radialAngleFromMList.Add(M[iM].M_ID, radialAngle);
                    }
                }
                M[cM].distanceFromMList = M[cM].distanceFromMList.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
            }

            List<int> inputM = M.Select(m => m.M_ID).ToList();

            int mPerTriplet = 3;

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> T = new List<mTriplet>();

            Image<Gray, Single> ofImg = new Image<Gray, Single>(1, 1);
            BlocksToPixels(ref ofImgBlock, ref ofImg, ofBlockSize);

            if (inputM.Count >= mPerTriplet)
            {
                for (int i = 0; i < inputM.Count; i++)
                {
                    int midx = Array.FindIndex(M, m => m.M_ID == inputM[i]);
                    var allNeighborID = M[midx].distanceFromMList.Keys.ToList();

                    /////Curve Axis, P=plus, M=Minus direction
                    List<PointF> SSampleP = new List<PointF>();
                    List<PointF> SSampleM = new List<PointF>();
                    List<Point> SSampleLineP = new List<Point>();
                    List<Point> SSampleLineM = new List<Point>();
                    List<Point> SSampleLineInMDirection = new List<Point>();
                    List<Point> SSampleLineInMDirectionOp = new List<Point>();
                    int l = 5;
                    int N = (int)(L / l) * 2;//number of iteration

                    OFFCForMinutia(ofImg, M[midx], L, N, false, SSampleM, SSampleLineM, l);
                    OFFCForMinutia(ofImg, M[midx], L, N, true, SSampleP, SSampleLineP, l);

                    if (SSampleM.Count == 1 && SSampleP.Count == 1)
                    {
                        double d = M[midx].direction;
                        double mD2OF = double.NaN;
                        if (d > Math.PI / 2.0 && d <= (3.0 * Math.PI / 2))
                        {
                            mD2OF = d - Math.PI;
                        }
                        else if (d > (3.0 * Math.PI / 2))
                        {
                            mD2OF = d - 2.0 * Math.PI;
                        }
                        else
                        {
                            mD2OF = d;
                        }

                        double djMinus = -1;
                        double djPlus = 1;
                        PointF start = new PointF(M[midx].p.X, M[midx].p.Y);
                        PointF endPlus = new PointF();
                        PointF endMinus = new PointF();
                        endPlus.X = (float)(start.X + (djPlus * l * Math.Cos(mD2OF)));
                        endPlus.Y = (float)(start.Y + (djPlus * l * Math.Sin(mD2OF)));
                        endMinus.X = (float)(start.X + (djMinus * l * Math.Cos(mD2OF)));
                        endMinus.Y = (float)(start.Y + (djMinus * l * Math.Sin(mD2OF)));

                        SSampleP.Add(new PointF(endPlus.X, endPlus.Y));
                        SSampleLineP.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endPlus)));
                        SSampleM.Add(new PointF(endMinus.X, endMinus.Y));
                        SSampleLineM.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endMinus)));
                    }

                    SSampleLineM = SSampleLineM.Distinct().ToList();
                    SSampleLineP = SSampleLineP.Distinct().ToList();

                    if (SSampleM.Count >= 2 || SSampleP.Count >= 2)
                    {
                        double ofM = ofImg.Data[(int)M[midx].p.Y, (int)M[midx].p.X, 0];
                        double difMDandOF = adPI(M[midx].direction, ofM);
                        double difMDandOF2 = adPI(ofM, M[midx].direction);

                        double polarAngle = Double.NaN;
                        bool MorP = true; // true = M
                        if (SSampleM.Count >= 2)
                        {
                            polarAngle = angleOf2Points(SSampleM[1], M[midx].p);
                            MorP = true;
                        }
                        else
                        {
                            polarAngle = angleOf2Points(SSampleP[1], M[midx].p);
                            MorP = false;
                        }
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, M[midx].direction);

                        if (Math.Abs(radialAngle) < Math.PI / 2)
                        {
                            if (MorP)
                            {
                                SSampleLineInMDirection = SSampleLineM;
                                SSampleLineInMDirectionOp = SSampleLineP;
                            }
                            else
                            {
                                SSampleLineInMDirection = SSampleLineP;
                                SSampleLineInMDirectionOp = SSampleLineM;
                            }
                        }
                        else
                        {
                            if (MorP)
                            {
                                SSampleLineInMDirection = SSampleLineP;
                                SSampleLineInMDirectionOp = SSampleLineM;
                            }
                            else
                            {
                                SSampleLineInMDirection = SSampleLineM;
                                SSampleLineInMDirectionOp = SSampleLineP;
                            }
                        }

                        Dictionary<int, List<double>> distanceToAxisInMDirection = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceOnAxisInMDirection = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceToAxisInMDirectionOp = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceOnAxisInMDirectionOp = new Dictionary<int, List<double>>();

                        Dictionary<int, double> neighborInMDirection = new Dictionary<int, double>();
                       
                        foreach (int mID in allNeighborID)
                        {
                            Point pInMDirection = new Point();
                            Point pInMDirectionOp = new Point();

                            PointF origin = new PointF(0, 0);
                            List<double> ARDistanceInMDirection = new List<double>();
                            List<double> ARDistanceInMDirectionOp = new List<double>();
                            List<PointF> wPOn = new List<PointF>();
                            List<PointF> wPOp = new List<PointF>();

                            foreach (Point p in SSampleLineInMDirection)
                            {
                                int mIDindex = Array.FindIndex(M, m => m.M_ID == mID);
                                double d = distanceOf2Points(M[mIDindex].p, p);
                                pInMDirection.X = p.X;
                                pInMDirection.Y = p.Y;
                                int idxP = SSampleLineInMDirection.IndexOf(pInMDirection);

                                if (distanceToAxisInMDirection.ContainsKey(mID))
                                {
                                    distanceToAxisInMDirection[mID].Add(d);
                                    distanceOnAxisInMDirection[mID].Add(idxP);
                                }
                                else
                                {
                                    List<double> tmpDTo = new List<double>();
                                    tmpDTo.Add(d);
                                    distanceToAxisInMDirection.Add(mID, tmpDTo);

                                    List<double> tmpDOn = new List<double>();
                                    tmpDOn.Add(idxP);
                                    distanceOnAxisInMDirection.Add(mID, tmpDOn);
                                }

                                PointF wP = new PointF((float)idxP, (float)d);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                ARDistanceInMDirection.Add(wScore);
                                wPOn.Add(wP);
                            }

                            foreach (Point p in SSampleLineInMDirectionOp)
                            {
                                int mIDindex = Array.FindIndex(M, m => m.M_ID == mID);
                                double d = distanceOf2Points(M[mIDindex].p, p);
                                pInMDirectionOp.X = p.X;
                                pInMDirectionOp.Y = p.Y;
                                int idxP = SSampleLineInMDirectionOp.IndexOf(pInMDirectionOp);

                                if (distanceToAxisInMDirectionOp.ContainsKey(mID))
                                {
                                    distanceToAxisInMDirectionOp[mID].Add(d);
                                    distanceOnAxisInMDirectionOp[mID].Add(idxP);
                                }
                                else
                                {
                                    List<double> tmpDTo = new List<double>();
                                    tmpDTo.Add(d);
                                    distanceToAxisInMDirectionOp.Add(mID, tmpDTo);

                                    List<double> tmpDOn = new List<double>();
                                    tmpDOn.Add(idxP);
                                    distanceOnAxisInMDirectionOp.Add(mID, tmpDOn);
                                }

                                PointF wP = new PointF((float)idxP, (float)d);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                ARDistanceInMDirectionOp.Add(wScore);
                                wPOp.Add(wP);
                            }

                            if (ARDistanceInMDirection.Min() <= ARDistanceInMDirectionOp.Min())
                            {
                                neighborInMDirection.Add(mID, ARDistanceInMDirection.Min());
                            }
                            else
                            {
                                int idxP = ARDistanceInMDirectionOp.IndexOf(ARDistanceInMDirectionOp.Min());
                                PointF wP = new PointF(-wPOp[idxP].X, wPOp[idxP].Y);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                neighborInMDirection.Add(mID, wScore);
                            }

                            wPOn.Clear();
                            wPOn.TrimExcess();
                            wPOn = null;

                            wPOp.Clear();
                            wPOp.TrimExcess();
                            wPOp = null;

                            ARDistanceInMDirectionOp.Clear();
                            ARDistanceInMDirectionOp.TrimExcess();
                            ARDistanceInMDirectionOp = null;

                            ARDistanceInMDirection.Clear();
                            ARDistanceInMDirection.TrimExcess();
                            ARDistanceInMDirection = null;
                        }

                        if (neighborInMDirection.Count >= n)
                        {
                            var id = neighborInMDirection.Keys.ToArray();
                            var score = neighborInMDirection.Values.ToArray();
                            Array.Sort(score, id);
                            var d = id.Take(n).ToList();

                            d.Insert(0, inputM[i]);

                            if (arrangeCWandremoveDuplicate)
                            {
                                #region arrangeClockwise M3gl
                                //M3gl
                                //List<int[]> tmp_mTriplets = new List<int[]>();
                                //tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(pID.ToArray(), mPerTriplet)); //same result as below for-loop
                                int ii = 0;
                                for (int j = 1; j < d.Count; j++)
                                {
                                    for (int k = j + 1; k < d.Count; k++)
                                    {
                                        int[] inmID = new int[3];
                                        inmID[0] = d[ii];
                                        inmID[1] = d[j];
                                        inmID[2] = d[k];

                                        int[] outmID = arrangeClockwise(M, inmID);

                                        pID_mTriplet.Add(outmID);
                                    }
                                }
                                #endregion
                            }
                            else
                            {
                                #region order the first neighboring minutia and the second neighboring minutia based on Jea and Govindaraju (2005)
                                List<int[]> tmp_mTriplets = new List<int[]>();
                                tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(d.ToArray(), mPerTriplet));
                                for (int mt = 0; mt < tmp_mTriplets.Count; mt++)
                                {
                                    int[] pID = new int[3];
                                    int idxP0 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][0]);
                                    int idxP1 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][1]);
                                    int idxP2 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][2]);

                                    //AR
                                    double signAcuteAngle = crossProductRightHandRule(M[idxP0].p, M[idxP1].p, M[idxP2].p);
                                    if (signAcuteAngle >= 0)
                                    {
                                        pID[0] = M[idxP0].M_ID;
                                        pID[1] = M[idxP1].M_ID;
                                        pID[2] = M[idxP2].M_ID;
                                    }
                                    else
                                    {
                                        pID[0] = M[idxP0].M_ID;
                                        pID[1] = M[idxP2].M_ID;
                                        pID[2] = M[idxP1].M_ID;
                                    }

                                    pID_mTriplet.Add(pID);
                                }
                                tmp_mTriplets.Clear();
                                tmp_mTriplets.TrimExcess();
                                tmp_mTriplets = null;
                                #endregion
                            }


                            d.Clear();
                            d.TrimExcess();
                            d = null;
                        }

                        neighborInMDirection.Clear();
                        neighborInMDirection = null;

                        distanceToAxisInMDirection.Clear();
                        distanceToAxisInMDirection = null;
                        distanceOnAxisInMDirection.Clear();
                        distanceOnAxisInMDirection = null;
                        distanceToAxisInMDirectionOp.Clear();
                        distanceToAxisInMDirectionOp = null;
                        distanceOnAxisInMDirectionOp.Clear();
                        distanceOnAxisInMDirectionOp = null;
                    }

                    SSampleP.Clear();
                    SSampleP.TrimExcess();
                    SSampleP = null;
                    SSampleM.Clear();
                    SSampleM.TrimExcess();
                    SSampleM = null;
                    SSampleLineP.Clear();
                    SSampleLineP.TrimExcess();
                    SSampleLineP = null;
                    SSampleLineM.Clear();
                    SSampleLineM.TrimExcess();
                    SSampleLineM = null;
                    SSampleLineInMDirection.Clear();
                    SSampleLineInMDirection.TrimExcess();
                    SSampleLineInMDirection = null;
                    SSampleLineInMDirectionOp.Clear();
                    SSampleLineInMDirectionOp.TrimExcess();
                    SSampleLineInMDirectionOp = null;

                    allNeighborID.Clear();
                    allNeighborID.TrimExcess();
                    allNeighborID = null;
                }

            }

            if (arrangeCWandremoveDuplicate)
            {
                pID_mTriplet = pID_mTriplet.Distinct(new IntArrayComparerIndexSensitive()).ToList();
            }

            ofImg.Dispose();

            for (int i = 0; i < pID_mTriplet.Count; i++)
            {
                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][0]);
                int idxP1 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][1]);
                int idxP2 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][2]);

                t.pi[0] = M[idxP0];
                t.pi[1] = M[idxP1];
                t.pi[2] = M[idxP2];

                T.Add(t);
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            inputM.Clear();
            inputM.TrimExcess();
            inputM = null;

            return T;
        }

        public static List<mTriplet> directionallyWeightedDistanceNearestNeighbor_Method1_34bw(KMinutia[] M, double D, double Theta, int n = 4, bool arrangeCWandremoveDuplicate = false)
        {
            for (int cM = 0; cM < M.Length; cM++)
            {
                M[cM].distanceFromMList = new Dictionary<int, double>();
                M[cM].radialAngleFromMList = new Dictionary<int, double>();
                M[cM].opRadialAngleFromMList = new Dictionary<int, double>();

                for (int iM = 0; iM < M.Length; iM++)
                {
                    if (cM != iM)
                    {
                        PointF p1 = M[cM].p;
                        double directionP1 = M[cM].direction;
                        double opDirectionP1 = M[cM].direction + Math.PI < 2 * Math.PI ? M[cM].direction + Math.PI : M[cM].direction - Math.PI;

                        PointF p2 = M[iM].p;
                        double directionP2 = M[iM].direction;

                        double distance = distanceOf2Points(p1, p2);
                        double polarAngle = angleOf2Points(p2, p1);
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, directionP1);
                        double opRadialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, opDirectionP1);
                        M[cM].distanceFromMList.Add(M[iM].M_ID, distance);
                        M[cM].radialAngleFromMList.Add(M[iM].M_ID, radialAngle);
                        M[cM].opRadialAngleFromMList.Add(M[iM].M_ID, opRadialAngle);
                    }
                }
                M[cM].distanceFromMList = M[cM].distanceFromMList.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
            }

            List<int> inputM = M.Select(m => m.M_ID).ToList();

            int mPerTriplet = 3;

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> T = new List<mTriplet>();
            if (inputM.Count >= mPerTriplet)
            {
                for (int i = 0; i < inputM.Count; i++)
                {
                    int midx = Array.FindIndex(M, m => m.M_ID == inputM[i]);

                    Dictionary<int, double> distance = new Dictionary<int, double>();
                    distance = M[inputM[i]].distanceFromMList.ToDictionary(k => k.Key, v => v.Value);
                    distance = distance.OrderBy(u => u.Value).ToDictionary(k => k.Key, v => v.Value);
                    var key = distance.Keys.ToArray();

                    var id = distance.Keys.ToArray();
                    var distanceList = distance.Values.ToArray();

                    double[] wRadius = new double[id.Length];
                    double[] wTheta = new double[id.Length];
                    double[] wScore = new double[id.Length];

                    for (int k = 0; k < id.Length; k++)
                    {
                        wRadius[k] = (distanceList[k] / D);
                        wTheta[k] = M[midx].radialAngleFromMList[id[k]] / Theta;
                        wScore[k] = Math.Sqrt(Math.Pow(wRadius[k], 2) + Math.Pow(wTheta[k], 2));
                    }
                    Array.Sort(wScore, id);

                    var d = id.Take(n - 2).ToList();

                    var idOp = distance.Keys.ToArray();
                    var distanceListOp = distance.Values.ToArray();
                    double[] wRadiusOp = new double[idOp.Length];
                    double[] wThetaOp = new double[idOp.Length];
                    double[] wScoreOp = new double[idOp.Length];

                    for (int k = 0; k < idOp.Length; k++)
                    {
                        wRadiusOp[k] = (distanceListOp[k] / D);
                        wThetaOp[k] = M[midx].opRadialAngleFromMList[idOp[k]] / Theta;
                        wScoreOp[k] = Math.Sqrt(Math.Pow(wRadiusOp[k], 2) + Math.Pow(wThetaOp[k], 2));
                    }
                    Array.Sort(wScoreOp, idOp);

                    for (int id34 = 0; id34 < idOp.Length; id34++)
                    {
                        if (d.Count < n)
                        {
                            if (!d.Contains(idOp[id34]))
                            {
                                d.Add(idOp[id34]);
                            }
                        }
                        else
                        {
                            break;
                        }
                    }

                    d.Insert(0, inputM[i]);

                    if (arrangeCWandremoveDuplicate)
                    {
                        #region arrangeClockwise M3gl
                        //M3gl
                        //List<int[]> tmp_mTriplets = new List<int[]>();
                        //tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(pID.ToArray(), mPerTriplet)); //same result as below for-loop
                        int ii = 0;
                        for (int j = 1; j < d.Count; j++)
                        {
                            for (int k = j + 1; k < d.Count; k++)
                            {
                                int[] inmID = new int[3];
                                inmID[0] = d[ii];
                                inmID[1] = d[j];
                                inmID[2] = d[k];

                                int[] outmID = arrangeClockwise(M, inmID);

                                pID_mTriplet.Add(outmID);
                            }
                        }
                        #endregion
                    }
                    else
                    {
                        #region order the first neighboring minutia and the second neighboring minutia based on Jea and Govindaraju (2005)
                        List<int[]> tmp_mTriplets = new List<int[]>();
                        tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(d.ToArray(), mPerTriplet));
                        for (int mt = 0; mt < tmp_mTriplets.Count; mt++)
                        {
                            int[] pID = new int[3];
                            int idxP0 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][0]);
                            int idxP1 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][1]);
                            int idxP2 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][2]);

                            //AR
                            double signAcuteAngle = crossProductRightHandRule(M[idxP0].p, M[idxP1].p, M[idxP2].p);
                            if (signAcuteAngle >= 0)
                            {
                                pID[0] = M[idxP0].M_ID;
                                pID[1] = M[idxP1].M_ID;
                                pID[2] = M[idxP2].M_ID;
                            }
                            else
                            {
                                pID[0] = M[idxP0].M_ID;
                                pID[1] = M[idxP2].M_ID;
                                pID[2] = M[idxP1].M_ID;
                            }

                            pID_mTriplet.Add(pID);
                        }
                        tmp_mTriplets.Clear();
                        tmp_mTriplets.TrimExcess();
                        tmp_mTriplets = null;
                        #endregion
                    }

                    d.Clear();
                    d.TrimExcess();
                    d = null;

                    distance.Clear();
                    distance = null;
                }
            }

            if (arrangeCWandremoveDuplicate)
            {
                pID_mTriplet = pID_mTriplet.Distinct(new IntArrayComparerIndexSensitive()).ToList();
            }

            for (int i = 0; i < pID_mTriplet.Count; i++)
            {
                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][0]);
                int idxP1 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][1]);
                int idxP2 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][2]);

                t.pi[0] = M[idxP0];
                t.pi[1] = M[idxP1];
                t.pi[2] = M[idxP2];

                T.Add(t);
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            inputM.Clear();
            inputM.TrimExcess();
            inputM = null;

            return T;
        }

        public static List<mTriplet> ridgeFlowDirectionallyWeightedDistanceNearestNeighbor_Method2_34bw(KMinutia[] M, Image<Gray, Single> ofImgBlock, int ofBlockSize, double D, double Theta, double L, int n = 4, bool arrangeCWandremoveDuplicate = false)
        {
            for (int cM = 0; cM < M.Length; cM++)
            {
                M[cM].distanceFromMList = new Dictionary<int, double>();
                M[cM].radialAngleFromMList = new Dictionary<int, double>();

                for (int iM = 0; iM < M.Length; iM++)
                {
                    if (cM != iM)
                    {
                        PointF p1 = M[cM].p;
                        double directionP1 = M[cM].direction;
                        double opDirectionP1 = M[cM].direction + Math.PI < 2 * Math.PI ? M[cM].direction + Math.PI : M[cM].direction - Math.PI;

                        PointF p2 = M[iM].p;
                        double directionP2 = M[iM].direction;

                        double distance = distanceOf2Points(p1, p2);
                        double polarAngle = angleOf2Points(p2, p1);
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, directionP1);
                        M[cM].distanceFromMList.Add(M[iM].M_ID, distance);
                        M[cM].radialAngleFromMList.Add(M[iM].M_ID, radialAngle);
                    }
                }
                M[cM].distanceFromMList = M[cM].distanceFromMList.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
            }

            List<int> inputM = M.Select(m => m.M_ID).ToList();

            int mPerTriplet = 3;

            List<int[]> pID_mTriplet = new List<int[]>();
            List<mTriplet> T = new List<mTriplet>();

            Image<Gray, Single> ofImg = new Image<Gray, Single>(1, 1);
            BlocksToPixels(ref ofImgBlock, ref ofImg, ofBlockSize);

            if (inputM.Count >= mPerTriplet)
            {
                for (int i = 0; i < inputM.Count; i++)
                {
                    int midx = Array.FindIndex(M, m => m.M_ID == inputM[i]);
                    var inRadius = M[midx].distanceFromMList.Keys.ToList();

                    /////Curve Axis
                    List<PointF> SSampleP = new List<PointF>();
                    List<PointF> SSampleM = new List<PointF>();
                    List<Point> SSampleLineP = new List<Point>();
                    List<Point> SSampleLineM = new List<Point>();
                    List<Point> SSampleLineInMDirection = new List<Point>();
                    List<Point> SSampleLineInMDirectionOp = new List<Point>();
                    int l = 5;
                    int N = (int)(L / l) * 2;

                    OFFCForMinutia(ofImg, M[midx], L, N, false, SSampleM, SSampleLineM, l);
                    OFFCForMinutia(ofImg, M[midx], L, N, true, SSampleP, SSampleLineP, l);

                    if (SSampleM.Count == 1 && SSampleP.Count == 1)
                    {
                        double d = M[midx].direction;
                        double mD2OF = double.NaN;
                        if (d > Math.PI / 2.0 && d <= (3.0 * Math.PI / 2))
                        {
                            mD2OF = d - Math.PI;
                        }
                        else if (d > (3.0 * Math.PI / 2))
                        {
                            mD2OF = d - 2.0 * Math.PI;
                        }
                        else
                        {
                            mD2OF = d;
                        }

                        double djMinus = -1;
                        double djPlus = 1;
                        PointF start = new PointF(M[midx].p.X, M[midx].p.Y);
                        PointF endPlus = new PointF();
                        PointF endMinus = new PointF();
                        endPlus.X = (float)(start.X + (djPlus * l * Math.Cos(mD2OF)));
                        endPlus.Y = (float)(start.Y + (djPlus * l * Math.Sin(mD2OF)));
                        endMinus.X = (float)(start.X + (djMinus * l * Math.Cos(mD2OF)));
                        endMinus.Y = (float)(start.Y + (djMinus * l * Math.Sin(mD2OF)));

                        SSampleP.Add(new PointF(endPlus.X, endPlus.Y));
                        SSampleLineP.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endPlus)));
                        SSampleM.Add(new PointF(endMinus.X, endMinus.Y));
                        SSampleLineM.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endMinus)));
                    }

                    SSampleLineM = SSampleLineM.Distinct().ToList();
                    SSampleLineP = SSampleLineP.Distinct().ToList();

                    if (SSampleM.Count >= 2 || SSampleP.Count >= 2)
                    {
                        double ofM = ofImg.Data[(int)M[midx].p.Y, (int)M[midx].p.X, 0];
                        double difMDandOF = adPI(M[midx].direction, ofM);
                        double difMDandOF2 = adPI(ofM, M[midx].direction);

                        double polarAngle = Double.NaN;
                        bool MorP = true; // true = M
                        if (SSampleM.Count >= 2)
                        {
                            polarAngle = angleOf2Points(SSampleM[1], M[midx].p);
                            MorP = true;
                        }
                        else
                        {
                            polarAngle = angleOf2Points(SSampleP[1], M[midx].p);
                            MorP = false;
                        }
                        double radialAngle = differentAngleFrom2Direction_0To2Pi(polarAngle, M[midx].direction);

                        if (Math.Abs(radialAngle) < Math.PI / 2)
                        {
                            if (MorP)
                            {
                                SSampleLineInMDirection = SSampleLineM;
                                SSampleLineInMDirectionOp = SSampleLineP;
                            }
                            else
                            {
                                SSampleLineInMDirection = SSampleLineP;
                                SSampleLineInMDirectionOp = SSampleLineM;
                            }
                        }
                        else
                        {
                            if (MorP)
                            {
                                SSampleLineInMDirection = SSampleLineP;
                                SSampleLineInMDirectionOp = SSampleLineM;
                            }
                            else
                            {
                                SSampleLineInMDirection = SSampleLineM;
                                SSampleLineInMDirectionOp = SSampleLineP;
                            }
                        }

                        Dictionary<int, List<double>> distanceToAxisInMDirection = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceOnAxisInMDirection = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceToAxisInMDirectionOp = new Dictionary<int, List<double>>();
                        Dictionary<int, List<double>> distanceOnAxisInMDirectionOp = new Dictionary<int, List<double>>();

                        Dictionary<int, double> neighborInMDirection = new Dictionary<int, double>();
                        Dictionary<int, double> neighborInMDirectionOp = new Dictionary<int, double>();

                        foreach (int mID in inRadius)
                        {
                            Point pInMDirection = new Point();
                            Point pInMDirectionOp = new Point();

                            PointF origin = new PointF(0, 0);
                            List<double> ARDistanceInMDirection = new List<double>();
                            List<double> ARDistanceInMDirectionOp = new List<double>();
                            List<PointF> wPOn = new List<PointF>();
                            List<PointF> wPOp = new List<PointF>();

                            foreach (Point p in SSampleLineInMDirection)
                            {
                                int mIDindex = Array.FindIndex(M, m => m.M_ID == mID);
                                double d = distanceOf2Points(M[mIDindex].p, p);
                                pInMDirection.X = p.X;
                                pInMDirection.Y = p.Y;
                                int idxP = SSampleLineInMDirection.IndexOf(pInMDirection);

                                if (distanceToAxisInMDirection.ContainsKey(mID))
                                {
                                    distanceToAxisInMDirection[mID].Add(d);
                                    distanceOnAxisInMDirection[mID].Add(idxP);
                                }
                                else
                                {
                                    List<double> tmpDTo = new List<double>();
                                    tmpDTo.Add(d);
                                    distanceToAxisInMDirection.Add(mID, tmpDTo);

                                    List<double> tmpDOn = new List<double>();
                                    tmpDOn.Add(idxP);
                                    distanceOnAxisInMDirection.Add(mID, tmpDOn);
                                }

                                PointF wP = new PointF((float)idxP, (float)d);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                ARDistanceInMDirection.Add(wScore);
                                wPOn.Add(wP);
                            }

                            foreach (Point p in SSampleLineInMDirectionOp)
                            {
                                int mIDindex = Array.FindIndex(M, m => m.M_ID == mID);
                                double d = distanceOf2Points(M[mIDindex].p, p);
                                pInMDirectionOp.X = p.X;
                                pInMDirectionOp.Y = p.Y;
                                int idxP = SSampleLineInMDirectionOp.IndexOf(pInMDirectionOp);

                                if (distanceToAxisInMDirectionOp.ContainsKey(mID))
                                {
                                    distanceToAxisInMDirectionOp[mID].Add(d);
                                    distanceOnAxisInMDirectionOp[mID].Add(idxP);
                                }
                                else
                                {
                                    List<double> tmpDTo = new List<double>();
                                    tmpDTo.Add(d);
                                    distanceToAxisInMDirectionOp.Add(mID, tmpDTo);

                                    List<double> tmpDOn = new List<double>();
                                    tmpDOn.Add(idxP);
                                    distanceOnAxisInMDirectionOp.Add(mID, tmpDOn);
                                }

                                PointF wP = new PointF((float)idxP, (float)d);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                ARDistanceInMDirectionOp.Add(wScore);
                                wPOp.Add(wP);
                            }

                            if (ARDistanceInMDirection.Min() <= ARDistanceInMDirectionOp.Min())
                            {
                                neighborInMDirection.Add(mID, ARDistanceInMDirection.Min());

                                //Op
                                int idxP = ARDistanceInMDirection.IndexOf(ARDistanceInMDirection.Min());
                                PointF wP = new PointF(-wPOn[idxP].X, wPOn[idxP].Y);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                neighborInMDirectionOp.Add(mID, wScore);
                            }
                            else
                            {
                                neighborInMDirectionOp.Add(mID, ARDistanceInMDirectionOp.Min());

                                int idxP = ARDistanceInMDirectionOp.IndexOf(ARDistanceInMDirectionOp.Min());
                                PointF wP = new PointF(-wPOp[idxP].X, wPOp[idxP].Y);
                                double distance = distanceOf2Points(origin, wP);
                                double polarAngleWP = angleOf2Points(wP, origin);
                                double radialAngleWP = differentAngleFrom2Direction_0To2Pi(polarAngleWP, 0);
                                double wRadius = distance / D;
                                double wTheta = radialAngleWP / Theta;
                                double wScore = Math.Sqrt(Math.Pow(wRadius, 2) + Math.Pow(wTheta, 2));

                                neighborInMDirection.Add(mID, wScore);
                            }

                            wPOn.Clear();
                            wPOn.TrimExcess();
                            wPOn = null;

                            wPOp.Clear();
                            wPOp.TrimExcess();
                            wPOp = null;

                            ARDistanceInMDirectionOp.Clear();
                            ARDistanceInMDirectionOp.TrimExcess();
                            ARDistanceInMDirectionOp = null;

                            ARDistanceInMDirection.Clear();
                            ARDistanceInMDirection.TrimExcess();
                            ARDistanceInMDirection = null;
                        }

                        List<int> selectedID = new List<int>();
                        bool formedTriplet = false;
                        if (neighborInMDirection.Count >= n - 2)
                        {
                            formedTriplet = true;
                            var id = neighborInMDirection.Keys.ToArray();
                            var score = neighborInMDirection.Values.ToArray();
                            Array.Sort(score, id);
                            selectedID.AddRange(id.Take(n - 2).ToList());
                        }

                        if (neighborInMDirectionOp.Count >= n - 2)
                        {
                            formedTriplet = true;
                            var idOp = neighborInMDirectionOp.Keys.ToArray();
                            var scoreOp = neighborInMDirectionOp.Values.ToArray();
                            Array.Sort(scoreOp, idOp);
                            for (int id34 = 0; id34 < idOp.Length; id34++)
                            {
                                if (selectedID.Count < n)
                                {
                                    if (!selectedID.Contains(idOp[id34]))
                                    {
                                        selectedID.Add(idOp[id34]);
                                    }
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }

                        if (formedTriplet)
                        {
                            selectedID.Insert(0, inputM[i]);

                            if (arrangeCWandremoveDuplicate)
                            {
                                #region arrangeClockwise M3gl
                                //M3gl
                                //List<int[]> tmp_mTriplets = new List<int[]>();
                                //tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(pID.ToArray(), mPerTriplet)); //same result as below for-loop
                                int ii = 0;
                                for (int j = 1; j < selectedID.Count; j++)
                                {
                                    for (int k = j + 1; k < selectedID.Count; k++)
                                    {
                                        int[] inmID = new int[3];
                                        inmID[0] = selectedID[ii];
                                        inmID[1] = selectedID[j];
                                        inmID[2] = selectedID[k];

                                        int[] outmID = arrangeClockwise(M, inmID);

                                        pID_mTriplet.Add(outmID);
                                    }
                                }
                                #endregion
                            }
                            else
                            {
                                #region order the first neighboring minutia and the second neighboring minutia based on Jea and Govindaraju (2005)
                                List<int[]> tmp_mTriplets = new List<int[]>();
                                tmp_mTriplets.AddRange(CombinationsRosettaWoRecursion(selectedID.ToArray(), mPerTriplet));
                                for (int mt = 0; mt < tmp_mTriplets.Count; mt++)
                                {
                                    int[] pID = new int[3];
                                    int idxP0 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][0]);
                                    int idxP1 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][1]);
                                    int idxP2 = Array.FindIndex(M, p => p.M_ID == tmp_mTriplets[mt][2]);

                                    //AR
                                    double signAcuteAngle = crossProductRightHandRule(M[idxP0].p, M[idxP1].p, M[idxP2].p);
                                    if (signAcuteAngle >= 0)
                                    {
                                        pID[0] = M[idxP0].M_ID;
                                        pID[1] = M[idxP1].M_ID;
                                        pID[2] = M[idxP2].M_ID;
                                    }
                                    else
                                    {
                                        pID[0] = M[idxP0].M_ID;
                                        pID[1] = M[idxP2].M_ID;
                                        pID[2] = M[idxP1].M_ID;
                                    }

                                    pID_mTriplet.Add(pID);
                                }
                                tmp_mTriplets.Clear();
                                tmp_mTriplets.TrimExcess();
                                tmp_mTriplets = null;
                                #endregion
                            }
                        }
                        selectedID.Clear();
                        selectedID.TrimExcess();
                        selectedID = null;

                        neighborInMDirection.Clear();
                        neighborInMDirection = null;
                        neighborInMDirectionOp.Clear();
                        neighborInMDirectionOp = null;

                        distanceToAxisInMDirection.Clear();
                        distanceToAxisInMDirection = null;
                        distanceOnAxisInMDirection.Clear();
                        distanceOnAxisInMDirection = null;
                        distanceToAxisInMDirectionOp.Clear();
                        distanceToAxisInMDirectionOp = null;
                        distanceOnAxisInMDirectionOp.Clear();
                        distanceOnAxisInMDirectionOp = null;
                    }

                    SSampleP.Clear();
                    SSampleP.TrimExcess();
                    SSampleP = null;
                    SSampleM.Clear();
                    SSampleM.TrimExcess();
                    SSampleM = null;
                    SSampleLineP.Clear();
                    SSampleLineP.TrimExcess();
                    SSampleLineP = null;
                    SSampleLineM.Clear();
                    SSampleLineM.TrimExcess();
                    SSampleLineM = null;
                    SSampleLineInMDirection.Clear();
                    SSampleLineInMDirection.TrimExcess();
                    SSampleLineInMDirection = null;
                    SSampleLineInMDirectionOp.Clear();
                    SSampleLineInMDirectionOp.TrimExcess();
                    SSampleLineInMDirectionOp = null;

                    inRadius.Clear();
                    inRadius.TrimExcess();
                    inRadius = null;
                }
            }

            if (arrangeCWandremoveDuplicate)
            {
                pID_mTriplet = pID_mTriplet.Distinct(new IntArrayComparerIndexSensitive()).ToList();
            }

            ofImg.Dispose();

            for (int i = 0; i < pID_mTriplet.Count; i++)
            {
                mTriplet t = new mTriplet();
                int idxP0 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][0]);
                int idxP1 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][1]);
                int idxP2 = Array.FindIndex(M, p => p.M_ID == pID_mTriplet[i][2]);

                t.pi[0] = M[idxP0];
                t.pi[1] = M[idxP1];
                t.pi[2] = M[idxP2];

                T.Add(t);
            }

            pID_mTriplet.Clear();
            pID_mTriplet.TrimExcess();
            pID_mTriplet = null;

            inputM.Clear();
            inputM.TrimExcess();
            inputM = null;

            return T;
        }

        public static void calDistanceOfPiToAllP(KMinutia[] P)
        {
            for (int pi = 0; pi < P.Length; pi++)
            {
                if (P[pi].distanceFromMList != null)
                {
                    P[pi].distanceFromMList.Clear();
                    P[pi].distanceFromMList = null;
                }
                P[pi].distanceFromMList = new Dictionary<int, double>();
                for (int pn = 0; pn < P.Length; pn++)
                {
                    if (pi != pn)
                    {
                        double distance = distanceOf2Points(P[pi].p, P[pn].p);
                        P[pi].distanceFromMList.Add(P[pn].M_ID, distance);
                    }
                }
                P[pi].distanceFromMList = P[pi].distanceFromMList.OrderBy(u => u.Value).ToDictionary(z => z.Key, y => y.Value);
            }

        }

        public static double distanceOf2Points(PointF a, PointF b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public static double angleOf2Points(PointF a, PointF b)
        {
            return Math.Atan2(a.Y - b.Y, a.X - b.X);
        }

        /// <summary>
        /// Calculate diffence between two angles.
        /// For angle measured increasing counter-clockwise, call this function as differentAngleFrom2Direction_0To2Pi(theta1, theta2).
        /// For angle measured increasing clockwise, call this function as differentAngleFrom2Direction_0To2Pi(theta2, theta1).
        /// The result is in range [-Pi, Pi].
        /// </summary>
        /// <param name="theta1">input angle1 in range 0 to 2Pi</param>
        /// <param name="theta2">input angle2 in range 0 to 2Pi</param>
        /// <returns>
        /// Return diffence of theta1 - theta2 if in range [-Pi, Pi] 
        /// or return 2Pi + diffence if the diffence less than -Pi,
        /// or return -2Pi + diffence if the diffence greater than Pi.
        /// </returns>
        public static double differentAngleFrom2Direction_0To2Pi(double theta1, double theta2)
        {
            double result = theta1 - theta2;
            if (result < -Math.PI)
            {
                result = (2.0 * Math.PI) + result;
            }
            else if (result >= Math.PI)
            {
                result = (-2.0 * Math.PI) + result;
            }
            return result;
        }

        public static int[] arrangeClockwise(KMinutia[] P, int[] mID)
        {
            int[] o = (int[])mID.Clone();
            double[] locationAngle = new double[mID.Length];
            int[] idxP = new int[mID.Length];

            idxP[0] = Array.FindIndex(P, p => p.M_ID == mID[0]);
            float CenterX = P[idxP[0]].p.X;
            float CenterY = P[idxP[0]].p.Y;
            for (int i = 1; i < mID.Length; i++)
            {
                idxP[i] = Array.FindIndex(P, p => p.M_ID == mID[i]);

                CenterX = (CenterX + P[idxP[i]].p.X) / 2;
                CenterY = (CenterY + P[idxP[i]].p.Y) / 2;
            }

            int idxP0 = Array.FindIndex(P, p => p.M_ID == mID[0]);
            int idxP1 = Array.FindIndex(P, p => p.M_ID == mID[1]);
            int idxP2 = Array.FindIndex(P, p => p.M_ID == mID[2]);

            float CenterX2 = (float)(((1.0 * P[idxP0].p.X + P[idxP1].p.X) / 2.0 + P[idxP2].p.X) / 2.0);
            float CenterY2 = (float)(((1.0 * P[idxP0].p.Y + P[idxP1].p.Y) / 2.0 + P[idxP2].p.Y) / 2.0);
            PointF centerP = new PointF(CenterX, CenterY);

            for (int i = 0; i < mID.Length; i++)
            {
                locationAngle[i] = ang(P[idxP[i]].p, centerP);
            }
            Array.Sort(locationAngle, o);

            return o;
        }

        // Enumerate all possible m-size combinations of [0, 1, ..., n-1] array
        // in lexicographic order (first [0, 1, 2, ..., m-1]).
        private static IEnumerable<int[]> CombinationsRosettaWoRecursion(int k, int n)
        {
            int[] result = new int[k];
            Stack<int> stack = new Stack<int>(k);
            stack.Push(0);
            while (stack.Count > 0)
            {
                int index = stack.Count - 1;
                int value = stack.Pop();
                while (value < n)
                {
                    result[index++] = value++;
                    stack.Push(value);
                    if (index != k) continue;
                    yield return (int[])result.Clone(); // thanks to @xanatos
                                                        //yield return result;
                    break;
                }
            }
        }

        public static IEnumerable<T[]> CombinationsRosettaWoRecursion<T>(T[] array, int k)
        {
            if (array.Length < k)
                throw new ArgumentException("Array length can't be less than number of selected elements");
            if (k < 1)
                throw new ArgumentException("Number of selected elements can't be less than 1");
            T[] result = new T[k];
            foreach (int[] j in CombinationsRosettaWoRecursion(k, array.Length))
            {
                for (int i = 0; i < k; i++)
                {
                    result[i] = array[j[i]];
                }
                //yield return (T[])result.Clone();
                if (result.Contains(array[0]))
                {
                    yield return (T[])result.Clone();
                }
            }
        }

        private static double crossProductRightHandRule(PointF p1, PointF p2, PointF p3)
        {
            double vPC1X = (p2.X - p1.X); double vPC1Y = (p2.Y - p1.Y);
            double vPC2X = (p3.X - p1.X); double vPC2Y = (p3.Y - p1.Y);

            double signAcuteAngle = (vPC1X * vPC2Y) - (vPC2X * vPC1Y);

            return signAcuteAngle;
        }

        /// <summary>
        /// ang(pi, pj) computes the angle of the vector with initial point at pi and terminal point at pj
        /// </summary>
        /// <param name="pi"></param>
        /// <param name="pj"></param>
        /// <returns></returns>
        public static double ang(PointF pi, PointF pj)
        {
            double result = 0.0;

            double deltaX = pi.X - pj.X;
            double deltaY = pi.Y - pj.Y;

            if (deltaX > 0 && deltaY >= 0)
            {
                result = Math.Atan(deltaY / deltaX);
            }
            else if (deltaX > 0 && deltaY < 0)
            {
                result = Math.Atan(deltaY / deltaX) + (2 * Math.PI);
            }
            else if (deltaX < 0)
            {
                result = Math.Atan(deltaY / deltaX) + Math.PI;
            }
            else if (deltaX == 0 && deltaY > 0)
            {
                result = Math.PI / 2;
            }
            else if (deltaX == 0 && deltaY < 0)
            {
                result = (3 * Math.PI) / 2;
            }

            return result;
        }

        /// <summary>
        /// adπ(alpha, beta) computes the minimum angle required to superpose two vectors with the same origin and angles alpha and beta respectively
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <returns>min(|alpha − beta|, 2π − |alpha − beta|)</returns>
        public static double adPI(double alpha, double beta)
        {
            double absResult = Math.Abs(alpha - beta);
            double _2PI_AbsResult = (2 * Math.PI) - Math.Abs(alpha - beta);

            return absResult < _2PI_AbsResult ? absResult : _2PI_AbsResult;
        }

        public static void OFFCForMinutia(Image<Gray, Single> ofImg, KMinutia M, double L, int N, bool dplus, List<PointF> SSample, List<Point> Line, int l = 5)
        {
            List<PointF> S = new List<PointF>();
            List<Point> SLine = new List<Point>();

            double djMinus = -1;
            double djPlus = 1;
            S.Add(new PointF(M.p.X, M.p.Y));

            bool first = true;
            double lastOF = Double.NaN;
            bool plus = dplus;
            bool radiusReached = false;

            while (!radiusReached && S.Count < N)
            {
                int x = (int)S[S.Count - 1].X;
                int y = (int)S[S.Count - 1].Y;

                double thetaS = ofImg.Data[y, x, 0];
                if (!Double.IsNaN(thetaS))
                {
                    PointF start = new PointF(S[S.Count - 1].X, S[S.Count - 1].Y);
                    PointF endPlus = new PointF();
                    PointF endMinus = new PointF();
                    endPlus.X = (float)(start.X + (djPlus * l * Math.Cos(thetaS)));
                    endPlus.Y = (float)(start.Y + (djPlus * l * Math.Sin(thetaS)));
                    endMinus.X = (float)(start.X + (djMinus * l * Math.Cos(thetaS)));
                    endMinus.Y = (float)(start.Y + (djMinus * l * Math.Sin(thetaS)));
                    LineSegment2DF sPlus = new LineSegment2DF(start, endPlus);
                    LineSegment2DF sMinus = new LineSegment2DF(start, endMinus);

                    if (!first)
                    {
                        double difOF = adPI(lastOF, thetaS);
                        plus = difOF >= ((90 * Math.PI) / 180) ? !plus : plus;
                        lastOF = thetaS;
                    }

                    if (plus)
                    {
                        float newXPlus = endPlus.X;
                        float newYPlus = endPlus.Y;

                        if (newXPlus < ofImg.Width && newYPlus < ofImg.Height && newXPlus >= 0 && newYPlus >= 0)
                        {
                            S.Add(new PointF(newXPlus, newYPlus));
                            SLine.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endPlus)));

                            radiusReached = distanceOf2Points(M.p, S[S.Count - 1]) < L ? false : true;
                        }
                        else
                        {
                            break;
                        }
                    }
                    else
                    {

                        float newXMinus = endMinus.X;
                        float newYMinus = endMinus.Y;

                        if (newXMinus < ofImg.Width && newYMinus < ofImg.Height && newXMinus >= 0 && newYMinus >= 0)
                        {
                            S.Add(new PointF(newXMinus, newYMinus));

                            SLine.AddRange(getBresenhamLinePoints(Point.Round(start), Point.Round(endMinus)));

                            radiusReached = distanceOf2Points(M.p, S[S.Count - 1]) < L ? false : true;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                else
                {
                    break;
                }

                if (first)
                {
                    first = false;
                    lastOF = thetaS;
                }
            }

            if (radiusReached)
            {
                for (int i = 0; i < SLine.Count; i++)
                {
                    if (distanceOf2Points(M.p, SLine[i]) > L)
                    {
                        SLine.RemoveAt(i);
                        i--;
                    }
                }
            }

            SSample.AddRange(S);
            Line.AddRange(SLine);

            S.Clear();
            S.TrimExcess();
            S = null;

            SLine.Clear();
            SLine.TrimExcess();
            SLine = null;
        }

        public static List<Point> getBresenhamLinePoints(Point p0, Point p1)
        {
            int x0 = p0.X;
            int y0 = p0.Y;
            int x1 = p1.X;
            int y1 = p1.Y;
            int dx = Math.Abs(x1 - x0);
            int dy = Math.Abs(y1 - y0);

            int sx = x0 < x1 ? 1 : -1;
            int sy = y0 < y1 ? 1 : -1;

            int err = dx - dy;

            var points = new List<Point>();

            while (true)
            {
                points.Add(new Point(x0, y0));
                if (x0 == x1 && y0 == y1) break;
                int e2 = 2 * err;
                if (e2 > -dy)
                {
                    err = err - dy;
                    x0 = x0 + sx;
                }
                if (e2 < dx)
                {
                    err = err + dx;
                    y0 = y0 + sy;
                }
            }

            return points;
        }

        public static void pixelToBlock(ref Image<Gray, Single> pixelImg, ref Image<Gray, Single> outBlockImg, int blockSize)
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
        }

        public static void BlocksToPixels(ref Image<Gray, Single> input, ref Image<Gray, Single> output, int blockSize)
        {
            Image<Gray, Single> pixelsImg = new Image<Gray, Single>(input.Width * blockSize, input.Height * blockSize);
            for (int inX = 0; inX < input.Width; inX++)
            {
                for (int inY = 0; inY < input.Height; inY++)
                {
                    float value = input.Data[inY, inX, 0];
                    for (int outX = inX * blockSize; outX < (inX * blockSize) + blockSize; outX++)
                    {
                        for (int outY = inY * blockSize; outY < (inY * blockSize) + blockSize; outY++)
                        {
                            pixelsImg.Data[outY, outX, 0] = value;
                        }
                    }
                }
            }

            output = pixelsImg.Clone();
            pixelsImg.Dispose();
        }

        public static void getOFFromFile(string ofTxtPath, Size fpImgSize, ref Image<Gray, Single> outOFImg, int ofBlockSize, Image<Gray, Single> fpSegmentImg = null)
        {
            string[] strOFLines = System.IO.File.ReadAllLines(ofTxtPath);
            string[] splitChars = { "\t", " ", "  " };

            int outputWidth = fpImgSize.Width / ofBlockSize;
            int outputHeight = fpImgSize.Height / ofBlockSize;
            if (fpImgSize.Width % ofBlockSize != 0)
            {
                outputWidth++;
            }
            if (fpImgSize.Height % ofBlockSize != 0)
            {
                outputHeight++;
            }

            Image<Gray, Single> ofImg = new Image<Gray, Single>(outputWidth, outputHeight);
            ofImg.SetZero();

            for (int y = 0; y < strOFLines.Length; y++)
            {
                string[] stringOF = strOFLines[y].Split(splitChars, StringSplitOptions.RemoveEmptyEntries);

                for (int x = 0; x < stringOF.Length; x++)
                {
                    if (fpSegmentImg == null)
                    {
                        ofImg.Data[y, x, 0] = (float)Convert.ToDouble(stringOF[x]);
                    }
                    else
                    {
                        if (fpSegmentImg.Data[y, x, 0] == 255)
                        {
                            ofImg.Data[y, x, 0] = (float)Convert.ToDouble(stringOF[x]);
                        }
                        else
                        {
                            ofImg[y, x] = new Gray(Double.NaN);
                        }
                    }
                }
            }

            outOFImg = ofImg.Clone();
            ofImg.Dispose();
        }

        public static void saveMTripletsToText(List<mTriplet> T, string savePath)
        {
            using (var writer = new StreamWriter(savePath))
            {
                for (int t = 0; t < T.Count; t++)
                {
                    writer.WriteLine("{0} {1} {2}", T[t].pi[0].M_ID, T[t].pi[1].M_ID, T[t].pi[2].M_ID);
                }
            }
        }

        public static void Quiver(ref Image<Gray, Single> Theta, ref Image<Gray, Single> Coh, ref Image<Gray, Single> Output, int blockSize, bool arrow = false)
        {
            double[] CohMin, CohMax;
            Point[] CohPMin, CohPMax;
            Coh.MinMax(out CohMin, out CohMax, out CohPMin, out CohPMax);

            Output = new Image<Gray, Single>(Theta.Width * blockSize, Theta.Height * blockSize);
            Output.SetValue(255);
            double maxLength = blockSize / 1.5;

            for (int y = 0; y < Theta.Height; y++)
            {
                for (int x = 0; x < Theta.Width; x++)
                {
                    double angle = Theta.Data[y, x, 0];
                    double length = (maxLength / CohMax[0]) * Coh.Data[y, x, 0];

                    if (!Double.IsNaN(angle))
                    {
                        PointF center = new PointF((x * blockSize) + (blockSize / 2), (y * blockSize) + (blockSize / 2));
                        PointF end = new PointF();
                        end.X = (float)(center.X + length * Math.Cos(angle));
                        end.Y = (float)(center.Y + length * Math.Sin(angle));

                        LineSegment2DF OFLine = new LineSegment2DF(center, end);
                        Output.Draw(OFLine, new Gray(0), 1);

                        if (arrow)
                        {
                            PointF pArrow1 = new PointF();
                            PointF pArrow2 = new PointF();
                            pArrow1.X = (float)Math.Round(((double)end.X - (length / 2) * Math.Cos(angle + Math.PI / 7)));
                            pArrow1.Y = (float)Math.Round(((double)end.Y - (length / 2) * Math.Sin(angle + Math.PI / 7)));
                            LineSegment2DF ArrowLine1 = new LineSegment2DF(end, pArrow1);
                            Output.Draw(ArrowLine1, new Gray(0), 1);
                            //QuiverBlock.Draw(ArrowLine1, new Gray(0), 1);

                            pArrow2.X = (float)Math.Round(((double)end.X - (length / 2) * Math.Cos(angle - Math.PI / 7)));
                            pArrow2.Y = (float)Math.Round(((double)end.Y - (length / 2) * Math.Sin(angle - Math.PI / 7)));
                            LineSegment2DF ArrowLine2 = new LineSegment2DF(end, pArrow2);
                            Output.Draw(ArrowLine2, new Gray(0), 1);
                            //QuiverBlock.Draw(ArrowLine2, new Gray(0), 1);
                        }
                    }
                    else
                    {
                        PointF center = new PointF((x * blockSize) + (blockSize / 2), (y * blockSize) + (blockSize / 2));
                        PointF end = center;

                        LineSegment2DF OFLine = new LineSegment2DF(center, end);
                        Output.Draw(OFLine, new Gray(0), 1);
                    }

                }
            }

        }

        public static void matchJYFNMR(string minutiaePath, string minutiaeTripletsPath, string fpImgPath, string fpImgExtension, string saveResultPath, string saveFilename, int numOfSamplePerFinger = 8, bool saveImgResult = false)
        {
            var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("*.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);
            
            string[] allFileNames = vAllFileNames.ToArray();

            Stopwatch stopwatch = new Stopwatch();
            TimeSpan totalElapsed = new TimeSpan();
            Stopwatch OneFingerStopwatch = new Stopwatch();

            Dictionary<string, MatchResult> finalScoreTable = new Dictionary<string, MatchResult>();

            //MTriplet_JY2000 jy = new MTriplet_JY2000();
            MTriplet_JY2000V2 jy = new MTriplet_JY2000V2();

            for (int g = 0; g < allFileNames.Length; g++)
            //for (int g = 0; g < allFileNames.Length - 2700; g++)//NIST14
            {
                OneFingerStopwatch.Restart();
                string mTripletPath = Path.GetDirectoryName(allFileNames[g]);
                string gfname = Path.GetFileNameWithoutExtension(allFileNames[g]);

                KMinutia[] J = readFingerNetMinutiae(allFileNames[g]);
                List<mTriplet> FlJ = readMinutiaeTripletsFromTextFile(ref J, minutiaeTripletsPath + gfname + ".txt");

                for (int q = g + 1; q < (g + numOfSamplePerFinger) - (g % numOfSamplePerFinger); q++)
                //for (int q = g + 2700; q < (g + 2700 + numOfSamplePerFinger - 1); q++)//NIST14
                {
                    double score = 0;

                    string qfname = Path.GetFileNameWithoutExtension(allFileNames[q]);

                    KMinutia[] I = readFingerNetMinutiae(allFileNames[q]);
                    List<mTriplet> FlI = readMinutiaeTripletsFromTextFile(ref I, minutiaeTripletsPath + qfname + ".txt");

                    List<int[]> MatchedPair = new List<int[]>();
                    stopwatch.Restart();
                    if (FlI.Count > 0 && FlJ.Count > 0)
                    {
                        score = jy.matchMinutiaeTriplets(FlI, FlJ, I, J, ref MatchedPair, false);
                    }
                    stopwatch.Stop();

                    totalElapsed = totalElapsed.Add(stopwatch.Elapsed);

                    MatchResult mr = new MatchResult();
                    mr.score = score;
                    mr.numOfMatch = MatchedPair.Count;
                    mr.elapsed = stopwatch.Elapsed;
                    mr.totalMilliseconds = stopwatch.Elapsed.TotalMilliseconds;
                    mr.ticks = stopwatch.Elapsed.Ticks;
                    finalScoreTable.Add(gfname + "VS" + qfname, mr);

                    #region plot matched minutiae
                    if (saveImgResult)
                    {
                        Image<Bgr, Single> galleryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + gfname + fpImgExtension);
                        Image<Bgr, Single> queryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + qfname + fpImgExtension);
                        Image<Bgr, Single> MatchedImage = plotMSTMatchedMinutiaeJY(galleryFPImage, queryFPImage, J, I, MatchedPair, FlJ, FlI);

                        if (!Directory.Exists(saveResultPath + saveFilename))
                        {
                            Directory.CreateDirectory(saveResultPath + saveFilename);
                        }

                        MatchedImage.Save(saveResultPath + saveFilename + "\\" + gfname + "VS" + qfname + "_" + MatchedPair.Count + ".png");
                        galleryFPImage.Dispose();
                        queryFPImage.Dispose();
                        MatchedImage.Dispose();
                    }
                    #endregion

                    MatchedPair.Clear();
                    MatchedPair.TrimExcess();
                    MatchedPair = null;

                    FlI.Clear();
                    FlI.TrimExcess();
                    FlI = null;
                }

                OneFingerStopwatch.Stop();
                Console.WriteLine(gfname + " has been matched against the remaining ones of the same finger. Time elapsed: {0}", OneFingerStopwatch.Elapsed);

                FlJ.Clear();
                FlJ.TrimExcess();
                FlJ = null;
            }
            saveFilename = saveFilename + ".txt";
            using (var writer = new StreamWriter(Path.Combine(saveResultPath, saveFilename)))
            {
                int r = 0;
                foreach (KeyValuePair<string, MatchResult> kv in finalScoreTable)
                {
                    writer.WriteLine("{0}\t{1}\t{2:R}\t{3}\t{4}\t{5}\t{6}", kv.Key, r += 1, kv.Value.score, kv.Value.numOfMatch, kv.Value.elapsed, kv.Value.totalMilliseconds, kv.Value.ticks);
                }
                //writer.WriteLine("Total time elapsed: {0}", totalElapsed);
                //writer.WriteLine("Total time elapsed (Milliseconds): {0}", totalElapsed.TotalMilliseconds);
                //writer.WriteLine("Total time elapsed (Ticks): {0}", totalElapsed.Ticks);
                //writer.WriteLine("Average time elapsed: {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison));
                //writer.WriteLine("Average time elapsed (Milliseconds): {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison).TotalMilliseconds);
            }

            finalScoreTable.Clear();
            finalScoreTable = null;
        }

        public static void matchJYFMR(string minutiaePath, string minutiaeTripletsPath, string fpImgPath, string fpImgExtension, string saveResultPath, string saveFilename, bool saveImgResult = false)
        {
            var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("*_1.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);
            //var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("S*.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);//NIST14

            string[] allFileNames = vAllFileNames.ToArray();

            Stopwatch stopwatch = new Stopwatch();
            TimeSpan totalElapsed = new TimeSpan();
            Stopwatch OneFingerStopwatch = new Stopwatch();

            Dictionary<string, MatchResult> finalScoreTable = new Dictionary<string, MatchResult>();

            //MTriplet_JY2000 jy = new MTriplet_JY2000();
            MTriplet_JY2000V2 jy = new MTriplet_JY2000V2();

            for (int g = 0; g < allFileNames.Length; g++)//107
            {
                OneFingerStopwatch.Restart();
                string mTripletPath = Path.GetDirectoryName(allFileNames[g]);
                string gfname = Path.GetFileNameWithoutExtension(allFileNames[g]);

                KMinutia[] J = readFingerNetMinutiae(allFileNames[g]);
                List<mTriplet> FlJ = readMinutiaeTripletsFromTextFile(ref J, minutiaeTripletsPath + gfname + ".txt");

                for (int q = g + 1; q < allFileNames.Length; q++)
                //for (int q = 0; q < allFileNames.Length; q++)//NIST14
                {
                    double score = 0;

                    string qfname = Path.GetFileNameWithoutExtension(allFileNames[q]);

                    //NIST14
                    //if (q == g)
                    //{
                    //    continue;
                    //}
                    //string qfPath = allFileNames[q].Replace(qfname, qfname.Replace("S", "F"));
                    //qfname = Path.GetFileNameWithoutExtension(qfPath);
                    //M3gl_Minutia[] I = readFingerNetMinutiae(qfPath);
                    //NIST14

                    KMinutia[] I = readFingerNetMinutiae(allFileNames[q]);
                    List<mTriplet> FlI = readMinutiaeTripletsFromTextFile(ref I, minutiaeTripletsPath + qfname + ".txt");

                    List<int[]> MatchedPair = new List<int[]>();
                    stopwatch.Restart();
                    if (FlI.Count > 0 && FlJ.Count > 0)
                    {
                        score = jy.matchMinutiaeTriplets(FlI, FlJ, I, J, ref MatchedPair, false);
                    }
                    stopwatch.Stop();

                    totalElapsed = totalElapsed.Add(stopwatch.Elapsed);

                    MatchResult mr = new MatchResult();
                    mr.score = score;
                    mr.numOfMatch = MatchedPair.Count;
                    mr.elapsed = stopwatch.Elapsed;
                    mr.totalMilliseconds = stopwatch.Elapsed.TotalMilliseconds;
                    mr.ticks = stopwatch.Elapsed.Ticks;
                    finalScoreTable.Add(qfname + "VS" + gfname, mr);

                    #region plot matched minutiae
                    if (saveImgResult)
                    {
                        Image<Bgr, Single> galleryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + gfname + fpImgExtension);
                        Image<Bgr, Single> queryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + qfname + fpImgExtension);
                        Image<Bgr, Single> MatchedImage = plotMSTMatchedMinutiaeJY(galleryFPImage, queryFPImage, J, I, MatchedPair, FlJ, FlI);

                        if (!Directory.Exists(saveResultPath + saveFilename))
                        {
                            Directory.CreateDirectory(saveResultPath + saveFilename);
                        }

                        MatchedImage.Save(saveResultPath + saveFilename + "\\" + gfname + "VS" + qfname + "_" + MatchedPair.Count + ".png");
                        galleryFPImage.Dispose();
                        queryFPImage.Dispose();
                        MatchedImage.Dispose();
                    }
                    #endregion

                    MatchedPair.Clear();
                    MatchedPair.TrimExcess();
                    MatchedPair = null;

                    FlI.Clear();
                    FlI.TrimExcess();
                    FlI = null;
                }

                OneFingerStopwatch.Stop();
                Console.WriteLine(gfname + " has been matched against the remaining ones of the same finger. Time elapsed: {0}", OneFingerStopwatch.Elapsed);

                FlJ.Clear();
                FlJ.TrimExcess();
                FlJ = null;
            }

            saveFilename = saveFilename + ".txt";
            using (var writer = new StreamWriter(Path.Combine(saveResultPath, saveFilename)))
            {
                int r = 0;
                foreach (KeyValuePair<string, MatchResult> kv in finalScoreTable)
                {
                    writer.WriteLine("{0}\t{1}\t{2:R}\t{3}\t{4}\t{5}\t{6}", kv.Key, r += 1, kv.Value.score, kv.Value.numOfMatch, kv.Value.elapsed, kv.Value.totalMilliseconds, kv.Value.ticks);
                }
                //writer.WriteLine("Total time elapsed: {0}", totalElapsed);
                //writer.WriteLine("Total time elapsed (Milliseconds): {0}", totalElapsed.TotalMilliseconds);
                //writer.WriteLine("Total time elapsed (Ticks): {0}", totalElapsed.Ticks);
                //writer.WriteLine("Average time elapsed: {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison));
                //writer.WriteLine("Average time elapsed (Milliseconds): {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison).TotalMilliseconds);
            }

            finalScoreTable.Clear();
            finalScoreTable = null;
        }

        public static void matchM3glFNMR(string minutiaePath, string minutiaeTripletsPath, string fpImgPath, string fpImgExtension, string saveResultPath, string saveFilename, int numOfSamplePerFinger = 8, bool saveImgResult = false)
        {
            var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("*.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);
            
            string[] allFileNames = vAllFileNames.ToArray();

            Stopwatch stopwatch = new Stopwatch();
            TimeSpan totalElapsed = new TimeSpan();
            Stopwatch OneFingerStopwatch = new Stopwatch();

            Dictionary<string, MatchResult> finalScoreTable = new Dictionary<string, MatchResult>();

            M3gl_Perez2012 m3gl = new M3gl_Perez2012();

            for (int g = 0; g < allFileNames.Length; g++)
            //for (int g = 0; g < allFileNames.Length - 2700; g++)//NIST14
            {
                OneFingerStopwatch.Restart();
                string mTripletPath = Path.GetDirectoryName(allFileNames[g]);
                string gfname = Path.GetFileNameWithoutExtension(allFileNames[g]);

                KMinutia[] P = readFingerNetMinutiae(allFileNames[g]);
                KMinutia[] allP = (KMinutia[])P.Clone();
                List<mTriplet> T = readMinutiaeTripletsFromTextFile(ref P, minutiaeTripletsPath + gfname + ".txt");


                for (int q = g + 1; q < (g + numOfSamplePerFinger) - (g % numOfSamplePerFinger); q++)//107
                //for (int q = g + 2700; q < (g + 2700 + numOfSamplePerFinger - 1); q++)//NIST14
                {
                    double score = 0;

                    string qfname = Path.GetFileNameWithoutExtension(allFileNames[q]);

                    KMinutia[] Q = readFingerNetMinutiae(allFileNames[q]);
                    KMinutia[] allQ = (KMinutia[])Q.Clone();
                    List<mTriplet> R = readMinutiaeTripletsFromTextFile(ref Q, minutiaeTripletsPath + qfname + ".txt");

                    List<int[]> MatchedPair = new List<int[]>();
                    stopwatch.Restart();

                    if (R.Count > 0 && T.Count > 0)
                    {
                        score = m3gl.matchMinutiaeTriplets(P, Q, T, R, ref MatchedPair);
                    }
                    stopwatch.Stop();

                    totalElapsed = totalElapsed.Add(stopwatch.Elapsed);

                    MatchResult mr = new MatchResult();
                    mr.score = score;
                    mr.numOfMatch = MatchedPair.Count;
                    mr.elapsed = stopwatch.Elapsed;
                    mr.totalMilliseconds = stopwatch.Elapsed.TotalMilliseconds;
                    mr.ticks = stopwatch.Elapsed.Ticks;
                    finalScoreTable.Add(gfname + "VS" + qfname, mr);

                    #region plot matched minutiae
                    if (saveImgResult)
                    {
                        Image<Bgr, Single> galleryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + gfname + fpImgExtension);
                        Image<Bgr, Single> queryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + qfname + fpImgExtension);
                        Image<Bgr, Single> MatchedImage = plotMSTMatchedMinutiaeM3gl(galleryFPImage, queryFPImage, allP, allQ, MatchedPair);

                        if (!Directory.Exists(saveResultPath + saveFilename))
                        {
                            Directory.CreateDirectory(saveResultPath + saveFilename);
                        }

                        MatchedImage.Save(saveResultPath + saveFilename + "\\" + gfname + "VS" + qfname + "_" + MatchedPair.Count + ".png");
                        galleryFPImage.Dispose();
                        queryFPImage.Dispose();
                        MatchedImage.Dispose();
                    }
                    #endregion

                    MatchedPair.Clear();
                    MatchedPair.TrimExcess();
                    MatchedPair = null;

                    R.Clear();
                    R.TrimExcess();
                    R = null;
                }

                OneFingerStopwatch.Stop();
                Console.WriteLine(gfname + " has been matched against the remaining ones of the same finger. Time elapsed: {0}", OneFingerStopwatch.Elapsed);

                T.Clear();
                T.TrimExcess();
                T = null;
            }

            saveFilename = saveFilename + ".txt";
            using (var writer = new StreamWriter(Path.Combine(saveResultPath, saveFilename)))
            {
                int r = 0;
                foreach (KeyValuePair<string, MatchResult> kv in finalScoreTable)
                {
                    writer.WriteLine("{0}\t{1}\t{2:R}\t{3}\t{4}\t{5}\t{6}", kv.Key, r += 1, kv.Value.score, kv.Value.numOfMatch, kv.Value.elapsed, kv.Value.totalMilliseconds, kv.Value.ticks);
                }
                //writer.WriteLine("Total time elapsed: {0}", totalElapsed);
                //writer.WriteLine("Total time elapsed (Milliseconds): {0}", totalElapsed.TotalMilliseconds);
                //writer.WriteLine("Total time elapsed (Ticks): {0}", totalElapsed.Ticks);
                //writer.WriteLine("Average time elapsed: {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison));
                //writer.WriteLine("Average time elapsed (Milliseconds): {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison).TotalMilliseconds);
            }

            finalScoreTable.Clear();
            finalScoreTable = null;
        }

        public static void matchM3glFMR(string minutiaePath, string minutiaeTripletsPath, string fpImgPath, string fpImgExtension, string saveResultPath, string saveFilename, bool saveImgResult = false)
        {
            var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("*_1.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);
            //var vAllFileNames = new DirectoryInfo(minutiaePath).GetFileSystemInfos("S*.mnt").OrderBy(fs => fs.Name, new NaturalStringComparer()).Select(x => x.FullName);//NIST14

            string[] allFileNames = vAllFileNames.ToArray();

            Stopwatch stopwatch = new Stopwatch();
            TimeSpan totalElapsed = new TimeSpan();
            Stopwatch OneFingerStopwatch = new Stopwatch();

            Dictionary<string, MatchResult> finalScoreTable = new Dictionary<string, MatchResult>();

            M3gl_Perez2012 m3gl = new M3gl_Perez2012();

            for (int g = 0; g < allFileNames.Length; g++)//107
            {
                OneFingerStopwatch.Restart();
                string mTripletPath = Path.GetDirectoryName(allFileNames[g]);
                string gfname = Path.GetFileNameWithoutExtension(allFileNames[g]);

                KMinutia[] P = readFingerNetMinutiae(allFileNames[g]);
                KMinutia[] allP = (KMinutia[])P.Clone();
                List<mTriplet> T = readMinutiaeTripletsFromTextFile(ref P, minutiaeTripletsPath + gfname + ".txt");


                for (int q = g + 1; q < allFileNames.Length; q++)//107
                //for (int q = 0; q < allFileNames.Length; q++)//NIST14
                {
                    double score = 0;

                    string qfname = Path.GetFileNameWithoutExtension(allFileNames[q]);

                    //NIST14
                    //if (q == g)
                    //{
                    //    continue;
                    //}
                    //string qfPath = allFileNames[q].Replace(qfname, qfname.Replace("S", "F"));
                    //qfname = Path.GetFileNameWithoutExtension(qfPath);
                    //M3gl_Minutia[] Q = readFingerNetMinutiae(qfPath);
                    //NIST14

                    KMinutia[] Q = readFingerNetMinutiae(allFileNames[q]);
                    KMinutia[] allQ = (KMinutia[])Q.Clone();
                    List<mTriplet> R = readMinutiaeTripletsFromTextFile(ref Q, minutiaeTripletsPath + qfname + ".txt");

                    List<int[]> MatchedPair = new List<int[]>();
                    stopwatch.Restart();

                    if (R.Count > 0 && T.Count > 0)
                    {
                        score = m3gl.matchMinutiaeTriplets(P, Q, T, R, ref MatchedPair);
                    }
                    stopwatch.Stop();

                    totalElapsed = totalElapsed.Add(stopwatch.Elapsed);

                    MatchResult mr = new MatchResult();
                    mr.score = score;
                    mr.numOfMatch = MatchedPair.Count;
                    mr.elapsed = stopwatch.Elapsed;
                    mr.totalMilliseconds = stopwatch.Elapsed.TotalMilliseconds;
                    mr.ticks = stopwatch.Elapsed.Ticks;
                    finalScoreTable.Add(qfname + "VS" + gfname, mr);

                    #region plot matched minutiae
                    if (saveImgResult)
                    {
                        Image<Bgr, Single> galleryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + gfname + fpImgExtension);
                        Image<Bgr, Single> queryFPImage = new Image<Bgr, Single>(fpImgPath + "\\" + qfname + fpImgExtension);
                        Image<Bgr, Single> MatchedImage = plotMSTMatchedMinutiaeM3gl(galleryFPImage, queryFPImage, allP, allQ, MatchedPair);

                        if (!Directory.Exists(saveResultPath + saveFilename))
                        {
                            Directory.CreateDirectory(saveResultPath + saveFilename);
                        }

                        MatchedImage.Save(saveResultPath + saveFilename + "\\" + gfname + "VS" + qfname + "_" + MatchedPair.Count + ".png");
                        galleryFPImage.Dispose();
                        queryFPImage.Dispose();
                        MatchedImage.Dispose();
                    }
                    #endregion

                    MatchedPair.Clear();
                    MatchedPair.TrimExcess();
                    MatchedPair = null;

                    R.Clear();
                    R.TrimExcess();
                    R = null;
                }

                OneFingerStopwatch.Stop();
                Console.WriteLine(gfname + " has been matched against the remaining ones of the same finger. Time elapsed: {0}", OneFingerStopwatch.Elapsed);

                T.Clear();
                T.TrimExcess();
                T = null;
            }

            saveFilename = saveFilename + ".txt";
            using (var writer = new StreamWriter(Path.Combine(saveResultPath, saveFilename)))
            {
                int r = 0;
                foreach (KeyValuePair<string, MatchResult> kv in finalScoreTable)
                {
                    writer.WriteLine("{0}\t{1}\t{2:R}\t{3}\t{4}\t{5}\t{6}", kv.Key, r += 1, kv.Value.score, kv.Value.numOfMatch, kv.Value.elapsed, kv.Value.totalMilliseconds, kv.Value.ticks);
                }
                //writer.WriteLine("Total time elapsed: {0}", totalElapsed);
                //writer.WriteLine("Total time elapsed (Milliseconds): {0}", totalElapsed.TotalMilliseconds);
                //writer.WriteLine("Total time elapsed (Ticks): {0}", totalElapsed.Ticks);
                //writer.WriteLine("Average time elapsed: {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison));
                //writer.WriteLine("Average time elapsed (Milliseconds): {0}", new TimeSpan(totalElapsed.Ticks / totalNumberOfComparison).TotalMilliseconds);
            }

            finalScoreTable.Clear();
            finalScoreTable = null;
        }

        public static Image<Bgr, Single> plotMSTMatchedMinutiaeJY(Image<Bgr, Single> galleryFPImage, Image<Bgr, Single> queryFPImage, KMinutia[] P, KMinutia[] Q, List<int[]> matchedPair, List<mTriplet> FlJ, List<mTriplet> FlI)
        {
            plotMinutiaeOnImage(ref galleryFPImage, P);
            plotMinutiaeOnImage(ref queryFPImage, Q);
            //ImageViewer.Show(galleryFPImage.ConcateHorizontal(queryFPImage));

            //Plot reference minutiae pair
            if (matchedPair.Count > 0)
            {
                int mIID = matchedPair[0][0];
                int mJID = matchedPair[0][1];
                var flI = FlI.Where(f => f.pi[0].M_ID == mIID).FirstOrDefault();
                var flJ = FlJ.Where(f => f.pi[0].M_ID == mJID).FirstOrDefault();
                int p0 = Array.FindIndex(P, m => m.M_ID == flJ.pi[0].M_ID);
                int p1 = Array.FindIndex(P, m => m.M_ID == flJ.pi[1].M_ID);
                int p2 = Array.FindIndex(P, m => m.M_ID == flJ.pi[2].M_ID);
                int q0 = Array.FindIndex(Q, m => m.M_ID == flI.pi[0].M_ID);
                int q1 = Array.FindIndex(Q, m => m.M_ID == flI.pi[1].M_ID);
                int q2 = Array.FindIndex(Q, m => m.M_ID == flI.pi[2].M_ID);

                PointF str = new PointF(P[p0].p.X, P[p0].p.Y);
                PointF end = new PointF();
                double length = 20;
                end.X = (float)(str.X + length * Math.Cos(P[p0].direction));
                end.Y = (float)(str.Y + length * Math.Sin(P[p0].direction));
                LineSegment2DF minutiaLine = new LineSegment2DF(str, end);
                galleryFPImage.Draw(new CircleF(P[p0].p, 5), new Bgr(Color.Lime), 2);
                galleryFPImage.Draw(minutiaLine, new Bgr(Color.Lime), 2);
                double lengthArrow = 7;
                PointF endArrow1 = new PointF();
                PointF endArrow2 = new PointF();
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(P[p0].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p0].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(P[p0].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p0].direction - (Math.PI / 7)));
                LineSegment2DF arrowLine1 = new LineSegment2DF(end, endArrow1);
                LineSegment2DF arrowLine2 = new LineSegment2DF(end, endArrow2);
                galleryFPImage.Draw(arrowLine1, new Bgr(Color.Lime), 2);
                galleryFPImage.Draw(arrowLine2, new Bgr(Color.Lime), 2);


                str = new PointF(P[p1].p.X, P[p1].p.Y);
                end.X = (float)(str.X + length * Math.Cos(P[p1].direction));
                end.Y = (float)(str.Y + length * Math.Sin(P[p1].direction));
                minutiaLine = new LineSegment2DF(str, end);
                galleryFPImage.Draw(new CircleF(P[p1].p, 5), new Bgr(Color.Magenta), 2);
                galleryFPImage.Draw(minutiaLine, new Bgr(Color.Magenta), 2);
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(P[p1].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p1].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(P[p1].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p1].direction - (Math.PI / 7)));
                arrowLine1 = new LineSegment2DF(end, endArrow1);
                arrowLine2 = new LineSegment2DF(end, endArrow2);
                galleryFPImage.Draw(arrowLine1, new Bgr(Color.Magenta), 2);
                galleryFPImage.Draw(arrowLine2, new Bgr(Color.Magenta), 2);


                str = new PointF(P[p2].p.X, P[p2].p.Y);
                end.X = (float)(str.X + length * Math.Cos(P[p2].direction));
                end.Y = (float)(str.Y + length * Math.Sin(P[p2].direction));
                minutiaLine = new LineSegment2DF(str, end);
                galleryFPImage.Draw(new CircleF(P[p2].p, 5), new Bgr(Color.Magenta), 2);
                galleryFPImage.Draw(minutiaLine, new Bgr(Color.Magenta), 2);
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(P[p2].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p2].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(P[p2].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(P[p2].direction - (Math.PI / 7)));
                arrowLine1 = new LineSegment2DF(end, endArrow1);
                arrowLine2 = new LineSegment2DF(end, endArrow2);
                galleryFPImage.Draw(arrowLine1, new Bgr(Color.Magenta), 2);
                galleryFPImage.Draw(arrowLine2, new Bgr(Color.Magenta), 2);


                str = new PointF(Q[q0].p.X, Q[q0].p.Y);
                end.X = (float)(str.X + length * Math.Cos(Q[q0].direction));
                end.Y = (float)(str.Y + length * Math.Sin(Q[q0].direction));
                minutiaLine = new LineSegment2DF(str, end);
                queryFPImage.Draw(new CircleF(Q[q0].p, 5), new Bgr(Color.Lime), 2);
                queryFPImage.Draw(minutiaLine, new Bgr(Color.Lime), 2);
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(Q[q0].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q0].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(Q[q0].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q0].direction - (Math.PI / 7)));
                arrowLine1 = new LineSegment2DF(end, endArrow1);
                arrowLine2 = new LineSegment2DF(end, endArrow2);
                queryFPImage.Draw(arrowLine1, new Bgr(Color.Lime), 2);
                queryFPImage.Draw(arrowLine2, new Bgr(Color.Lime), 2);


                str = new PointF(Q[q1].p.X, Q[q1].p.Y);
                end.X = (float)(str.X + length * Math.Cos(Q[q1].direction));
                end.Y = (float)(str.Y + length * Math.Sin(Q[q1].direction));
                minutiaLine = new LineSegment2DF(str, end);
                queryFPImage.Draw(new CircleF(Q[q1].p, 5), new Bgr(Color.Magenta), 2);
                queryFPImage.Draw(minutiaLine, new Bgr(Color.Magenta), 2);
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(Q[q1].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q1].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(Q[q1].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q1].direction - (Math.PI / 7)));
                arrowLine1 = new LineSegment2DF(end, endArrow1);
                arrowLine2 = new LineSegment2DF(end, endArrow2);
                queryFPImage.Draw(arrowLine1, new Bgr(Color.Magenta), 2);
                queryFPImage.Draw(arrowLine2, new Bgr(Color.Magenta), 2);


                str = new PointF(Q[q2].p.X, Q[q2].p.Y);
                end.X = (float)(str.X + length * Math.Cos(Q[q2].direction));
                end.Y = (float)(str.Y + length * Math.Sin(Q[q2].direction));
                minutiaLine = new LineSegment2DF(str, end);
                queryFPImage.Draw(new CircleF(Q[q2].p, 5), new Bgr(Color.Magenta), 2);
                queryFPImage.Draw(minutiaLine, new Bgr(Color.Magenta), 2);
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(Q[q2].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q2].direction + (Math.PI / 7)));
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(Q[q2].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(Q[q2].direction - (Math.PI / 7)));
                arrowLine1 = new LineSegment2DF(end, endArrow1);
                arrowLine2 = new LineSegment2DF(end, endArrow2);
                queryFPImage.Draw(arrowLine1, new Bgr(Color.Magenta), 2);
                queryFPImage.Draw(arrowLine2, new Bgr(Color.Magenta), 2);
            }

            if (matchedPair.Count > 0)
            {
                Random r = new Random((int)(DateTime.Now.Ticks & 0x0000ffff));
                Bgr[] matchedLineColor = new Bgr[matchedPair.Count];
                int c = 0;
                while (true)
                {
                    Bgr color = new Bgr(r.NextDouble() * 255 + 90, r.NextDouble() * 255 + 90, r.NextDouble() * 255 + 90);
                    matchedLineColor[c] = color;
                    c++;

                    if (c == matchedLineColor.Length)
                    {
                        break;
                    }

                }

                int[] MP = new int[matchedPair.Count];
                int[] MQ = new int[matchedPair.Count];

                List<KMinutia> matchedMInP = new List<KMinutia>();

                for (int mp = 0; mp < matchedPair.Count; mp++)
                {
                    MP[mp] = matchedPair[mp][1];
                    MQ[mp] = matchedPair[mp][0];

                    int mpidx = Array.FindIndex(P, m => m.M_ID == matchedPair[mp][1]);
                    KMinutia mtmp = new KMinutia();
                    mtmp.M_ID = P[mpidx].M_ID;
                    mtmp.p = P[mpidx].p;
                    mtmp.direction = P[mpidx].direction;
                    matchedMInP.Add(mtmp);
                }

                MatchMinutiaeCluster cluster = new MatchMinutiaeCluster();
                cluster.m = matchedMInP.ToArray();
                calMSTFromClusterMinutiae(ref cluster);

                int eline = 0;
                foreach (Edge<int> e in cluster.mstPrim)
                {
                    if (e != null)
                    {
                        int matchedidx1 = Array.FindIndex(MP, v => v == e.From.Data);
                        int matchedidx2 = Array.FindIndex(MP, v => v == e.To.Data);

                        int mPidx1 = Array.FindIndex(P, v => v.M_ID == MP[matchedidx1]);
                        int mPidx2 = Array.FindIndex(P, v => v.M_ID == MP[matchedidx2]);
                        PointF mP1 = P[mPidx1].p;
                        PointF mP2 = P[mPidx2].p;

                        int mQidx1 = Array.FindIndex(Q, v => v.M_ID == MQ[matchedidx1]);
                        int mQidx2 = Array.FindIndex(Q, v => v.M_ID == MQ[matchedidx2]);
                        PointF mQ1 = Q[mQidx1].p;
                        PointF mQ2 = Q[mQidx2].p;

                        LineSegment2DF mPLine = new LineSegment2DF(mP1, mP2);
                        LineSegment2DF mQLine = new LineSegment2DF(mQ1, mQ2);
                        galleryFPImage.Draw(mPLine, matchedLineColor[eline], 2);
                        queryFPImage.Draw(mQLine, matchedLineColor[eline], 2);
                    }
                    eline++;
                }

                matchedMInP.Clear();
                matchedMInP.TrimExcess();
                matchedMInP = null;
            }

            Image<Bgr, Single> imgMatched = queryFPImage.ConcateHorizontal(galleryFPImage).Clone();
            //ImageViewer.Show(imgMatched);

            return imgMatched;
        }

        public static Image<Bgr, Single> plotMSTMatchedMinutiaeM3gl(Image<Bgr, Single> galleryFPImage, Image<Bgr, Single> queryFPImage, KMinutia[] P, KMinutia[] Q, List<int[]> matchedPair)
        {
            plotMinutiaeOnImage(ref galleryFPImage, P);
            plotMinutiaeOnImage(ref queryFPImage, Q);

            if (matchedPair.Count > 0)
            {
                Random r = new Random((int)(DateTime.Now.Ticks & 0x0000ffff));
                Bgr[] matchedLineColor = new Bgr[matchedPair.Count];
                int c = 0;
                while (true)
                {
                    Bgr color = new Bgr(r.NextDouble() * 255 + 90, r.NextDouble() * 255 + 90, r.NextDouble() * 255 + 90);
                    matchedLineColor[c] = color;
                    c++;

                    if (c == matchedLineColor.Length)
                    {
                        break;
                    }

                }

                int[] MP = new int[matchedPair.Count];
                int[] MQ = new int[matchedPair.Count];

                List<KMinutia> matchedMInP = new List<KMinutia>();

                for (int mp = 0; mp < matchedPair.Count; mp++)
                {
                    MP[mp] = matchedPair[mp][1];
                    MQ[mp] = matchedPair[mp][0];

                    int mpidx = Array.FindIndex(P, m => m.M_ID == matchedPair[mp][1]);
                    KMinutia mtmp = new KMinutia();
                    mtmp.M_ID = P[mpidx].M_ID;
                    mtmp.p = P[mpidx].p;
                    mtmp.direction = P[mpidx].direction;
                    matchedMInP.Add(mtmp);
                }

                MatchMinutiaeCluster cluster = new MatchMinutiaeCluster();
                cluster.m = matchedMInP.ToArray();
                calMSTFromClusterMinutiae(ref cluster);

                int eline = 0;
                foreach (Edge<int> e in cluster.mstPrim)
                {
                    if (e != null)
                    {
                        int matchedidx1 = Array.FindIndex(MP, v => v == e.From.Data);
                        int matchedidx2 = Array.FindIndex(MP, v => v == e.To.Data);

                        int mPidx1 = Array.FindIndex(P, v => v.M_ID == MP[matchedidx1]);
                        int mPidx2 = Array.FindIndex(P, v => v.M_ID == MP[matchedidx2]);
                        PointF mP1 = P[mPidx1].p;
                        PointF mP2 = P[mPidx2].p;

                        int mQidx1 = Array.FindIndex(Q, v => v.M_ID == MQ[matchedidx1]);
                        int mQidx2 = Array.FindIndex(Q, v => v.M_ID == MQ[matchedidx2]);
                        PointF mQ1 = Q[mQidx1].p;
                        PointF mQ2 = Q[mQidx2].p;

                        LineSegment2DF mPLine = new LineSegment2DF(mP1, mP2);
                        LineSegment2DF mQLine = new LineSegment2DF(mQ1, mQ2);
                        galleryFPImage.Draw(mPLine, matchedLineColor[eline], 2);
                        queryFPImage.Draw(mQLine, matchedLineColor[eline], 2);
                    }
                    eline++;
                }

                matchedMInP.Clear();
                matchedMInP.TrimExcess();
                matchedMInP = null;
            }

            Image<Bgr, Single> imgMatched = queryFPImage.ConcateHorizontal(galleryFPImage).Clone();
            //ImageViewer.Show(imgMatched);

            return imgMatched;
        }

        public static void plotMinutiaeOnImage(ref Image<Bgr, Single> FPImage, KMinutia[] M)
        {
            for (int m = 0; m < M.Length; m++)
            {
                PointF str = new PointF(M[m].p.X, M[m].p.Y);
                PointF end = new PointF();
                double length = 15;
                end.X = (float)(str.X + length * Math.Cos(M[m].direction));
                end.Y = (float)(str.Y + length * Math.Sin(M[m].direction));
                LineSegment2DF minutiaLine = new LineSegment2DF(str, end);
                FPImage.Draw(new CircleF(M[m].p, 5), new Bgr(0, 0, 255), 1);
                FPImage.Draw(minutiaLine, new Bgr(0, 0, 255), 1);
                //FPImage.Draw(M[m].M_ID.ToString(), new Point((int)str.X - 5, (int)str.Y + 13), FontFace.HersheyPlain, 0.5, new Bgr(0, 0, 255));

                length = 20;
                end.X = (float)(str.X + length * Math.Cos(M[m].direction));
                end.Y = (float)(str.Y + length * Math.Sin(M[m].direction));
                minutiaLine = new LineSegment2DF(str, end);
                FPImage.Draw(new CircleF(M[m].p, 5), new Bgr(0, 0, 255), 2);
                FPImage.Draw(minutiaLine, new Bgr(0, 0, 255), 2);

                double lengthArrow = 7;
                PointF endArrow1 = new PointF();
                endArrow1.X = (float)(end.X - lengthArrow * Math.Cos(M[m].direction + (Math.PI / 7)));
                endArrow1.Y = (float)(end.Y - lengthArrow * Math.Sin(M[m].direction + (Math.PI / 7)));
                PointF endArrow2 = new PointF();
                endArrow2.X = (float)(end.X - lengthArrow * Math.Cos(M[m].direction - (Math.PI / 7)));
                endArrow2.Y = (float)(end.Y - lengthArrow * Math.Sin(M[m].direction - (Math.PI / 7)));
                LineSegment2DF arrowLine1 = new LineSegment2DF(end, endArrow1);
                LineSegment2DF arrowLine2 = new LineSegment2DF(end, endArrow2);
                FPImage.Draw(arrowLine1, new Bgr(Color.Red), 2);
                FPImage.Draw(arrowLine2, new Bgr(Color.Red), 2);
            }
        }

        #region Minimum Spanning Tree
        public static void calMSTFromClusterMinutiae(ref MatchMinutiaeCluster cluster)
        {
            Triangle2DF[] delaunayTriangles;
            VoronoiFacet[] voronoiFacets;
            var clusterM = cluster.m.Select(id => id.M_ID).ToList();
            var ptsInCluster = cluster.m.Select(p => (PointF)p.p).ToArray();
            CreateSubdivision(ptsInCluster, out delaunayTriangles, out voronoiFacets);

            cluster.gCluster = new Graph<int>(false, true);
            List<Node<int>> node = new List<Node<int>>();

            cluster.gMSTPrim = new Graph<int>(false, false);
            cluster.gMSTPrim._isAttributed = true;
            List<Node<int>> nodeMSTPrim = new List<Node<int>>();
            foreach (int m in clusterM)
            {
                node.Add(cluster.gCluster.AddNode(m));

                nodeMSTPrim.Add(cluster.gMSTPrim.AddNode(m));
            }

            List<Node<int>> isolateTriangleNode = new List<Node<int>>();
            foreach (Triangle2DF dt in delaunayTriangles)
            {
                if (ptsInCluster.Contains(dt.V0) && ptsInCluster.Contains(dt.V1) && ptsInCluster.Contains(dt.V2))
                {
                    var mID0 = cluster.m.Where(m => m.p == dt.V0).Select(m => m.M_ID).FirstOrDefault();
                    var mID1 = cluster.m.Where(m => m.p == dt.V1).Select(m => m.M_ID).FirstOrDefault();
                    var mID2 = cluster.m.Where(m => m.p == dt.V2).Select(m => m.M_ID).FirstOrDefault();

                    int idx0 = clusterM.IndexOf(mID0);
                    int idx1 = clusterM.IndexOf(mID1);
                    int idx2 = clusterM.IndexOf(mID2);

                    int w0 = (int)(distanceOf2Points(cluster.m[idx0].p, cluster.m[idx1].p) + 0.5);
                    int w1 = (int)(distanceOf2Points(cluster.m[idx0].p, cluster.m[idx2].p) + 0.5);
                    int w2 = (int)(distanceOf2Points(cluster.m[idx1].p, cluster.m[idx2].p) + 0.5);

                    if (cluster.gCluster[idx0, idx1] == null)
                    {
                        cluster.gCluster.AddEdge(node[idx0], node[idx1], w0);
                    }
                    if (cluster.gCluster[idx0, idx2] == null)
                    {
                        cluster.gCluster.AddEdge(node[idx0], node[idx2], w1);
                    }
                    if (cluster.gCluster[idx1, idx2] == null)
                    {
                        cluster.gCluster.AddEdge(node[idx1], node[idx2], w2);
                    }

                    #region Check Isolate Triangle
                    var v0Count = delaunayTriangles.Where(d => d.V0 == dt.V0 || d.V1 == dt.V0 || d.V2 == dt.V0).Count();
                    var v1Count = delaunayTriangles.Where(d => d.V0 == dt.V1 || d.V1 == dt.V1 || d.V2 == dt.V1).Count();
                    var v2Count = delaunayTriangles.Where(d => d.V0 == dt.V2 || d.V1 == dt.V2 || d.V2 == dt.V2).Count();

                    if (v0Count == 1 && v1Count == 1 && v2Count == 1)
                    {
                        isolateTriangleNode.Add(node[idx0]);
                        isolateTriangleNode.Add(node[idx1]);
                        isolateTriangleNode.Add(node[idx2]);
                    }
                    #endregion

                }
            }

            addEdgeForUnconnectedNode(isolateTriangleNode, ref cluster);

            var unconnectedList = cluster.gCluster.Nodes.Where(n => n.Neighbors.Count == 0).ToList();
            addEdgeForUnconnectedNode(unconnectedList, ref cluster);

            if (cluster.gCluster.Nodes.Count > 0)
            {
                cluster.mstPrim = cluster.gCluster.MinimumSpanningTreePrim();
            }
            else
            {
                cluster.mstPrim = new List<Edge<int>>();
            }

            isolateTriangleNode.Clear();
            isolateTriangleNode.TrimExcess();
            isolateTriangleNode = null;
            unconnectedList.Clear();
            unconnectedList.TrimExcess();
            unconnectedList = null;
            node.Clear();
            node.TrimExcess();
            node = null;
            nodeMSTPrim.Clear();
            nodeMSTPrim.TrimExcess();
            nodeMSTPrim = null;
            clusterM.Clear();
            clusterM.TrimExcess();
            clusterM = null;
        }

        /// <summary>
        /// Create planar subdivision for random points
        /// </summary>
        public static void CreateSubdivision(PointF[] pts, out Triangle2DF[] delaunayTriangles, out VoronoiFacet[] voronoiFacets)
        {
            using (Subdiv2D subdivision = new Subdiv2D(pts))
            {
                //Obtain the delaunay's triangulation from the set of points;
                delaunayTriangles = subdivision.GetDelaunayTriangles();

                //Obtain the voronoi facets from the set of points
                voronoiFacets = subdivision.GetVoronoiFacets();
            }
        }

        public static void addEdgeForUnconnectedNode(List<Node<int>> UnconnectedNode, ref MatchMinutiaeCluster cluster)
        {
            var clusterM = cluster.gCluster.Nodes.Select(n => n.Data).ToList();
            foreach (Node<int> n in UnconnectedNode)
            {
                int mIDN = n.Data;
                int idxN = clusterM.IndexOf(mIDN);

                foreach (Node<int> c in cluster.gCluster.Nodes)
                {
                    int mIDC = c.Data;
                    int idxC = clusterM.IndexOf(mIDC);
                    if (n.Data != c.Data && cluster.gCluster[idxN, idxC] == null)
                    {
                        int wN = (int)(distanceOf2Points(cluster.m[idxN].p, cluster.m[idxC].p) + 0.5);
                        cluster.gCluster.AddEdge(cluster.gCluster.Nodes[idxN], cluster.gCluster.Nodes[idxC], wN);
                    }
                }
            }
            clusterM.Clear();
            clusterM.TrimExcess();
            clusterM = null;
        }
        #endregion
    }

    [SuppressUnmanagedCodeSecurity]
    internal static class SafeNativeMethods
    {
        [DllImport("shlwapi.dll", CharSet = CharSet.Unicode)]
        public static extern int StrCmpLogicalW(string psz1, string psz2);
    }

    public sealed class NaturalStringComparer : IComparer<string>
    {
        public int Compare(string a, string b)
        {
            return SafeNativeMethods.StrCmpLogicalW(a, b);
        }
    }

    class IntArrayComparerIndexSensitive : IEqualityComparer<int[]>
    {
        public bool Equals(int[] x, int[] y)
        {
            //Check whether the compared objects reference the same data.
            if (Object.ReferenceEquals(x, y)) return true;

            int numOfSameData = 0;
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] == y[i])
                {
                    numOfSameData++;
                }
            }
            if (numOfSameData == x.Length)
            {
                return true;
            }

            // If got this far, arrays are not equal
            return false;

        }

        public int GetHashCode(int[] intArray)
        {
            //Check whether the object is null
            if (Object.ReferenceEquals(intArray, null)) return 0;

            //Calculate the hash code for the array
            int hashCode = 0;
            bool isFirst = true;
            foreach (int i in intArray)
            {
                if (isFirst)
                {
                    hashCode = i;
                    isFirst = false;
                }
                else
                {
                    hashCode = hashCode ^ i;
                }
            }
            return hashCode;
        }
    }

    public class KMinutia
    {
        public int M_ID;
        public Point p;
        public double direction;
        public string type;

        public Dictionary<int, double> distanceFromMList;
        public Dictionary<int, double> radialAngleFromMList;
        public Dictionary<int, double> opRadialAngleFromMList;
        public bool valid;
    }

    public class mTriplet
    {
        public KMinutia[] pi;
        public double[] di;
        public double dMax;
        public double dMid;
        public double dMin;
        public double[] alphai;
        public double[] betai;
        public int[] ridgeCounts;

        public mTriplet()
        {
            pi = new KMinutia[3];
            di = new double[3];
            alphai = new double[6];
            betai = new double[3];
        }
    };

    public class Pose
    {
        public PointF p;
        public double angle;
    }

    public class MatchResult
    {
        public string ID;
        public int rank;
        public double score;
        public int numOfMatch;

        public TimeSpan elapsed;
        public double totalMilliseconds;
        public long ticks;
    }

    public class MatchMinutiaeCluster
    {
        public KMinutia[] m;
        public Graph<int> gCluster;
        public List<Edge<int>> mstPrim;
        public Graph<int> gMSTPrim;
    }
}
