using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//For PointF
using System.Drawing;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class MTriplet_JY2000
    {
        static int numberOfFeatures = 11;
        double bl = 6 * numberOfFeatures;
        double Wd = 1;
        double WThetaPhi = (0.3 * 180) / Math.PI;
        double Wn = 3;
        double Wt = 3;
        double BgD = 8;
        double BgTheta = Math.PI / 6;
        double BgPhi = Math.PI / 6;

        public void calIntrinsicFeatures(ref mTriplet r)
        {
            r.di = new double[2];
            r.di[0] = EuclideanDistanceOf2Points(r.pi[0].p, r.pi[1].p);
            r.di[1] = EuclideanDistanceOf2Points(r.pi[0].p, r.pi[2].p);

            r.alphai = new double[2];
            r.alphai[0] = Difference2Pi(r.pi[0].direction, ang(r.pi[0].p, r.pi[1].p));
            r.alphai[1] = Difference2Pi(r.pi[0].direction, ang(r.pi[0].p, r.pi[2].p));

            //double alphai0 = differentAngleFrom2Direction_0To2Pi(angleOf2Points(r.pi[1].p, r.pi[0].p), r.pi[0].direction);
            //double alphai1 = differentAngleFrom2Direction_0To2Pi(angleOf2Points(r.pi[2].p, r.pi[0].p), r.pi[0].direction);
            //r.alphai[0] = alphai0;
            //r.alphai[1] = alphai1;

            r.betai = new double[2];
            r.betai[0] = Difference2Pi(r.pi[0].direction, r.pi[1].direction);
            r.betai[1] = Difference2Pi(r.pi[0].direction, r.pi[2].direction);

            //double betai0 = differentAngleFrom2Direction_0To2Pi(r.pi[0].direction, r.pi[1].direction);
            //double betai1 = differentAngleFrom2Direction_0To2Pi(r.pi[0].direction, r.pi[2].direction);
            //r.betai[0] = betai0;
            //r.betai[1] = betai1;
        }

        public double matchMinutiaeTriplets(List<mTriplet> FlI, List<mTriplet> FlJ, KMinutia[] I, KMinutia[] J, ref List<int[]> MatchedPair, bool RidgeNType)
        {
            List<int[]> posMatchedPairs = new List<int[]>();
            Dictionary<string, double> sl = new Dictionary<string, double>();
            Dictionary<string, double> ml = new Dictionary<string, double>();

            List<int> matchedIid = new List<int>();
            List<int> matchedJid = new List<int>();
            Dictionary<int, double> matchedScore = new Dictionary<int, double>();
            List<int[]> finalMatchedPairs = new List<int[]>();

            for (int i = 0; i < FlI.Count; i++)
            {
                for (int j = 0; j < FlJ.Count; j++)
                {
                    double tmp_sl = matchLocalFeatures(FlI[i], FlJ[j], RidgeNType);

                    sl.Add(FlI[i].pi[0].M_ID.ToString() + "_" + FlJ[j].pi[0].M_ID.ToString(), tmp_sl);

                    int[] mPair = new int[2];
                    mPair[0] = FlI[i].pi[0].M_ID;
                    mPair[1] = FlJ[j].pi[0].M_ID;
                    posMatchedPairs.Add(mPair);
                }
            }

            sl = sl.OrderByDescending(s => s.Value).ToDictionary(k => k.Key, v => v.Value);

            string b1_b2 = sl.ElementAt(0).Key;
            string[] b1b2 = b1_b2.Split('_');
            int b1 = Convert.ToInt32(b1b2[0]);
            int b2 = Convert.ToInt32(b1b2[1]);
            int idxb1 = Array.FindIndex(I, p => p.M_ID == b1);
            int idxb2 = Array.FindIndex(J, p => p.M_ID == b2);

            matchedIid.Add(b1);
            matchedJid.Add(b2);
            double b1b2_ml = 0.5 + (0.5 * sl.ElementAt(0).Value);
            matchedScore.Add(b1, b1b2_ml);

            for (int m = 0; m < posMatchedPairs.Count; m++)
            {
                int i = posMatchedPairs[m][0];
                int j = posMatchedPairs[m][1];
                string keyMatchID = i.ToString() + "_" + j.ToString();

                int idxi = Array.FindIndex(I, p => p.M_ID == i);
                int idxj = Array.FindIndex(J, p => p.M_ID == j);

                double EDib1 = EuclideanDistanceOf2Points(I[idxi].p, I[idxb1].p);
                double EDjb2 = EuclideanDistanceOf2Points(J[idxj].p, J[idxb2].p);

                double matchedED = EDib1 - EDjb2;

                double alphaib1 = Difference2Pi(I[idxb1].direction, ang(I[idxi].p, I[idxb1].p));
                double alphajb2 = Difference2Pi(J[idxb2].direction, ang(J[idxj].p, J[idxb2].p));
                double diff = Math.Abs(alphaib1 - alphajb2);
                double matchedAlpha = Math.Min(diff, 2 * Math.PI - diff);

                //double alphaib1 = differentAngleFrom2Direction_0To2Pi(angleOf2Points(I[idxi].p, I[idxb1].p), I[idxb1].direction);
                //double alphajb2 = differentAngleFrom2Direction_0To2Pi(angleOf2Points(J[idxj].p, J[idxb2].p), J[idxb2].direction);
                //double matchedAlpha = differentAngleFrom2Direction_0To2Pi(alphaib1, alphajb2);

                double betaib1 = Difference2Pi(I[idxb1].direction, I[idxi].direction);
                double betajb2 = Difference2Pi(J[idxb2].direction, J[idxj].direction);
                diff = Math.Abs(betaib1 - betajb2);
                double matchedBeta = Math.Min(diff, 2 * Math.PI - diff);

                //double betaib1 = differentAngleFrom2Direction_0To2Pi(I[idxb1].direction, I[idxi].direction);
                //double betajb2 = differentAngleFrom2Direction_0To2Pi(J[idxb2].direction, J[idxj].direction);
                //double matchedBeta = differentAngleFrom2Direction_0To2Pi(betaib1, betajb2);

                if (Math.Abs(matchedED) < BgD && Math.Abs(matchedAlpha) < BgTheta && Math.Abs(matchedBeta) < BgPhi)
                {
                    double tmp_ml = 0.5 + (0.5 * sl[keyMatchID]);
                    if (matchedIid.Contains(i) || matchedJid.Contains(j))
                    {
                        int previousIdxi = matchedIid.FindIndex(mid => mid.Equals(i));
                        int previousIdxj = matchedJid.FindIndex(mid => mid.Equals(j));
                        if (previousIdxi != -1 && previousIdxj != -1)
                        {
                            if (previousIdxi == previousIdxj)
                            {
                                int keyI = matchedIid[previousIdxi];
                                double previous_ml = matchedScore[keyI];
                                if (tmp_ml > previous_ml)
                                {
                                    matchedIid.RemoveAt(previousIdxi);
                                    matchedJid.RemoveAt(previousIdxi);
                                    matchedScore.Remove(keyI);

                                    matchedIid.Add(i);
                                    matchedJid.Add(j);
                                    matchedScore.Add(i, tmp_ml);
                                }
                            }
                            else
                            {
                                int keyI1 = matchedIid[previousIdxi];
                                double previous_ml1 = matchedScore[keyI1];

                                int keyI2 = matchedIid[previousIdxj];
                                double previous_ml2 = matchedScore[keyI2];

                                if (tmp_ml > previous_ml1 && tmp_ml > previous_ml2)
                                {
                                    matchedIid.RemoveAt(previousIdxi);
                                    matchedJid.RemoveAt(previousIdxi);
                                    matchedScore.Remove(keyI1);

                                    int idxjAfterremoveIdxi = previousIdxj > previousIdxi ? previousIdxj - 1 : previousIdxj;
                                    matchedIid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedJid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedScore.Remove(keyI2);

                                    matchedIid.Add(i);
                                    matchedJid.Add(j);
                                    matchedScore.Add(i, tmp_ml);
                                }
                                else if (tmp_ml > previous_ml1)
                                {
                                    matchedIid.RemoveAt(previousIdxi);
                                    matchedJid.RemoveAt(previousIdxi);
                                    matchedScore.Remove(keyI1);

                                    int idxjAfterremoveIdxi = previousIdxj > previousIdxi ? previousIdxj - 1 : previousIdxj;
                                    matchedIid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedJid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedScore.Remove(keyI2);

                                    matchedIid.Add(i);
                                    matchedJid.Add(j);
                                    matchedScore.Add(i, tmp_ml);
                                }
                                else if (tmp_ml > previous_ml2)
                                {
                                    matchedIid.RemoveAt(previousIdxi);
                                    matchedJid.RemoveAt(previousIdxi);
                                    matchedScore.Remove(keyI1);

                                    int idxjAfterremoveIdxi = previousIdxj > previousIdxi ? previousIdxj - 1 : previousIdxj;
                                    matchedIid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedJid.RemoveAt(idxjAfterremoveIdxi);
                                    matchedScore.Remove(keyI2);

                                    matchedIid.Add(i);
                                    matchedJid.Add(j);
                                    matchedScore.Add(i, tmp_ml);
                                }
                            }
                        }
                        else if (previousIdxi != -1)
                        {
                            int keyI = matchedIid[previousIdxi];
                            double previous_ml = matchedScore[keyI];
                            if (tmp_ml > previous_ml)
                            {
                                matchedIid.RemoveAt(previousIdxi);
                                matchedJid.RemoveAt(previousIdxi);
                                matchedScore.Remove(keyI);

                                matchedIid.Add(i);
                                matchedJid.Add(j);
                                matchedScore.Add(i, tmp_ml);
                            }
                        }
                        else if (previousIdxj != -1)
                        {
                            int keyI = matchedIid[previousIdxj];
                            double previous_ml = matchedScore[keyI];
                            if (tmp_ml > previous_ml)
                            {
                                matchedIid.RemoveAt(previousIdxj);
                                matchedJid.RemoveAt(previousIdxj);
                                matchedScore.Remove(keyI);

                                matchedIid.Add(i);
                                matchedJid.Add(j);
                                matchedScore.Add(i, tmp_ml);
                            }
                        }
                    }
                    else
                    {
                        matchedIid.Add(i);
                        matchedJid.Add(j);
                        matchedScore.Add(i, tmp_ml);
                    }
                }

            }

            double score = 0;
            for (int i = 0; i < matchedIid.Count; i++)
            {
                int[] matchedPair = new int[2];
                matchedPair[0] = matchedIid[i];
                matchedPair[1] = matchedJid[i];
                finalMatchedPairs.Add(matchedPair);
                score += matchedScore[matchedIid[i]];
            }
            score = 100 * (score / Math.Max(I.Length, J.Length));

            MatchedPair = finalMatchedPairs;

            matchedScore.Clear();
            matchedScore = null;

            matchedJid.Clear();
            matchedJid.TrimExcess();
            matchedJid = null;
            matchedIid.Clear();
            matchedIid.TrimExcess();
            matchedIid = null;

            ml.Clear();
            ml = null;
            sl.Clear();
            sl = null;

            posMatchedPairs.Clear();
            posMatchedPairs.TrimExcess();
            posMatchedPairs = null;

            return score;
        }

        public double matchLocalFeatures(mTriplet fli, mTriplet flj, bool RidgeNType)
        {
            double score = 0;

            calIntrinsicFeatures(ref fli);
            calIntrinsicFeatures(ref flj);

            double di0 = fli.di[0] - flj.di[0];
            double di1 = fli.di[1] - flj.di[1];

            double alphai0 = differentAngleFrom2Direction_0To2Pi(fli.alphai[0], flj.alphai[0]);
            double alphai1 = differentAngleFrom2Direction_0To2Pi(fli.alphai[1], flj.alphai[1]);

            double betai0 = differentAngleFrom2Direction_0To2Pi(fli.betai[0], flj.betai[0]);
            double betai1 = differentAngleFrom2Direction_0To2Pi(fli.betai[1], flj.betai[1]);

            double sl = 0;

            if (RidgeNType)
            {
                double rc0 = fli.ridgeCounts[0] - flj.ridgeCounts[0];
                double rc1 = fli.ridgeCounts[1] - flj.ridgeCounts[1];

                double t0 = fli.pi[0].type == flj.pi[0].type ? 0 : 1;
                double t1 = fli.pi[1].type == flj.pi[1].type ? 0 : 1;
                double t2 = fli.pi[2].type == flj.pi[2].type ? 0 : 1;

                numberOfFeatures = 11;
                bl = 6 * numberOfFeatures;

                //double magnitudeDiff = Math.Sqrt((Wd * Math.Pow(di0, 2)) + (Wd * Math.Pow(di1, 2))
                //                                + (WThetaPhi * Math.Pow(alphai0, 2)) + (WThetaPhi * Math.Pow(alphai1, 2))
                //                                + (WThetaPhi * Math.Pow(betai0, 2)) + (WThetaPhi * Math.Pow(betai1, 2))
                //                                + (Wn * Math.Pow(rc0, 2)) + (Wn * Math.Pow(rc1, 2))
                //                                + (Wt * Math.Pow(t0, 2)) + (Wt * Math.Pow(t1, 2)) + (Wt * Math.Pow(t2, 2)));

                double absDiff = (Wd * Math.Abs(di0)) + (Wd * Math.Abs(di1))
                                                + (WThetaPhi * Math.Abs(alphai0)) + (WThetaPhi * Math.Abs(alphai1))
                                                + (WThetaPhi * Math.Abs(betai0)) + (WThetaPhi * Math.Abs(betai1))
                                                + (Wn * Math.Abs(rc0)) + (Wn * Math.Abs(rc1))
                                                + (Wt * Math.Abs(t0)) + (Wt * Math.Abs(t1)) + (Wt * Math.Abs(t2));

                sl = absDiff;//magnitudeDiff;
            }
            else
            {
                numberOfFeatures = 6;
                bl = 6 * numberOfFeatures;

                //double magnitudeDiff = Math.Sqrt((Wd * Math.Pow(di0, 2)) + (Wd * Math.Pow(di1, 2))
                //                                + (WThetaPhi * Math.Pow(alphai0, 2)) + (WThetaPhi * Math.Pow(alphai1, 2))
                //                                + (WThetaPhi * Math.Pow(betai0, 2)) + (WThetaPhi * Math.Pow(betai1, 2)));

                double absDiff = (Wd * Math.Abs(di0)) + (Wd * Math.Abs(di1))
                                                + (WThetaPhi * Math.Abs(alphai0)) + (WThetaPhi * Math.Abs(alphai1))
                                                + (WThetaPhi * Math.Abs(betai0)) + (WThetaPhi * Math.Abs(betai1));

                sl = absDiff;//magnitudeDiff;
            }

            score = sl < bl ? (bl - sl) / bl : 0;

            return score;
        }

        public double EuclideanDistanceOf2Points(Point a, Point b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public double Difference2Pi(double alpha, double beta)
        {
            if (beta >= alpha)
                return (beta - alpha);
            return beta - alpha + 2 * Math.PI;
        }
        
        public double ang(Point pi, Point pj)
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

        public double angleOf2Points(PointF a, PointF b)
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
        public double differentAngleFrom2Direction_0To2Pi(double theta1, double theta2)
        {
            double result = theta1 - theta2;
            if (result <= -Math.PI)
            {
                result = (2.0 * Math.PI) + result;
            }
            else if (result > Math.PI)
            {
                result = (2.0 * Math.PI) - result;
            }
            return result;
        }
    }
}
