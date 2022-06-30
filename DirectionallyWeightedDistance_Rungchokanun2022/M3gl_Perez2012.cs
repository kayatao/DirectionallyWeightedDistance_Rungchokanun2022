using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//For PointF
using System.Drawing;

namespace DirectionallyWeightedDistance_Rungchokanun2022
{
    class M3gl_Perez2012
    {
        double tl = 12;
        double tg = 12;
        double ta = Math.PI / 6;

        public double matchMinutiaeTriplets(KMinutia[] P, KMinutia[] Q, List<mTriplet> T, List<mTriplet> R, ref List<int[]> MatchedPair)
        {
            List<int[]> M = localMinutiaeMatching(T, R);
            int n = globalMinutiaeMatching(M, Q, P, ref MatchedPair);
            double score = Math.Pow(n, 2) / (P.Length * Q.Length);

            M.Clear();
            M.TrimExcess();
            M = null;

            return score;
        }

        public List<int[]> localMinutiaeMatching(List<mTriplet> T, List<mTriplet> R)
        {
            List<mTripletPair> A = new List<mTripletPair>();
            List<int[]> M = new List<int[]>();

            T.ForEach(t => calIntrinsicFeatures(ref t));
            R.ForEach(r => calIntrinsicFeatures(ref r));

            foreach (mTriplet ri in R)
            {
                var Tk = T.Where(t => Math.Abs(t.dMax - ri.dMax) <= tl && Math.Abs(t.dMid - ri.dMid) <= tl && Math.Abs(t.dMin - ri.dMin) <= tl).ToList();

                foreach (mTriplet ti in Tk)
                {
                    mTriplet rMax = new mTriplet();
                    double similarity = sInv(ti, ri, ref rMax);
                    if (similarity != 0)
                    {
                        mTripletPair mtp = new mTripletPair();
                        mtp.r = rMax;
                        mtp.t = ti;
                        mtp.similarity = similarity;
                        A.Add(mtp);
                    }
                }
            }

            A.Sort((a, b) => b.similarity.CompareTo(a.similarity));

            List<int> Mq = new List<int>();
            List<int> Mp = new List<int>();
            foreach (mTripletPair mtp in A)
            {
                int[][] B = new int[3][];
                B[0] = new int[2];
                B[1] = new int[2];
                B[2] = new int[2];

                B[0][0] = mtp.r.pi[0].M_ID;
                B[0][1] = mtp.t.pi[0].M_ID;
                B[1][0] = mtp.r.pi[1].M_ID;
                B[1][1] = mtp.t.pi[1].M_ID;
                B[2][0] = mtp.r.pi[2].M_ID;
                B[2][1] = mtp.t.pi[2].M_ID;

                foreach (int[] b in B)
                {
                    if (!(Mq.Contains(b[0]) || Mp.Contains(b[1])))
                    {
                        Mq.Add(b[0]);
                        Mp.Add(b[1]);
                        int[] m = new int[2];
                        m[0] = b[0];
                        m[1] = b[1];
                        M.Add(m);
                    }
                }
            }

            Mq.Clear();
            Mq.TrimExcess();
            Mq = null;
            Mp.Clear();
            Mp.TrimExcess();
            Mp = null;

            return M;

        }

        public int globalMinutiaeMatching(List<int[]> M, KMinutia[] Q, KMinutia[] P, ref List<int[]> MatchedPair)
        {
            int n = 0;
            int qRef = -1;
            int pRef = -1;

            for (int i = 0; i < M.Count; i++)
            {
                List<int> Eq = new List<int>();
                List<int> Ep = new List<int>();
                List<int[]> E = new List<int[]>();

                int qiID = M[i][0];
                int piID = M[i][1];

                int idx_qiID = Array.FindIndex(Q, q => q.M_ID == qiID);
                int idx_piID = Array.FindIndex(P, p => p.M_ID == piID);

                KMinutia qi = Q[idx_qiID];
                KMinutia pi = P[idx_piID];
                int x1 = qi.p.X;
                int y1 = qi.p.Y;
                double theta1 = qi.direction;
                int x2 = pi.p.X;
                int y2 = pi.p.Y;
                double theta2 = pi.direction;

                double deltaTheta = theta2 - theta1;
                double cosDeltaTheta = Math.Cos(deltaTheta);
                double sinDeltaTheta = Math.Sin(deltaTheta);

                Eq.Add(qiID);
                Ep.Add(piID);
                int[] ei = new int[2];
                ei[0] = qiID;
                ei[1] = piID;
                E.Add(ei);

                for (int j = 0; j < M.Count; j++)
                {
                    int qjID = M[j][0];
                    int pjID = M[j][1];
                    if (!Eq.Contains(qjID) && !Ep.Contains(pjID))
                    {
                        int idx_qjID = Array.FindIndex(Q, q => q.M_ID == qjID);
                        int idx_pjID = Array.FindIndex(P, p => p.M_ID == pjID);

                        KMinutia qj = Q[idx_qjID];
                        KMinutia pj = P[idx_pjID];
                        int x3 = qj.p.X;
                        int y3 = qj.p.Y;
                        double theta3 = qj.direction;
                        int x4 = pj.p.X;
                        int y4 = pj.p.Y;
                        double theta4 = pj.direction;

                        double x_ = (cosDeltaTheta * (x3 - x1)) + (-sinDeltaTheta * (y3 - y1)) + x2;
                        double y_ = (sinDeltaTheta * (x3 - x1)) + (cosDeltaTheta * (y3 - y1)) + y2;
                        double theta_ = (theta3 - theta1) + theta2;
                        theta_ = (theta_ > 2 * Math.PI) ? theta_ - 2 * Math.PI : (theta_ < 0) ? theta_ + 2 * Math.PI : theta_;

                        PointF q_ = new PointF((float)x_, (float)y_);
                        PointF p_ = new PointF((float)x4, (float)y4);
                        double alignedDistance = EuclideanDistanceOf2Points(q_, p_);
                        double theta21 = ad2PI(theta2, theta1);
                        double theta43 = ad2PI(theta4, theta3);
                        double difPairedTheta = adPI(theta21, theta43);
                        double difAlignedTheta = adPI(theta_, theta4);

                        if (alignedDistance <= tg && difPairedTheta <= ta && difAlignedTheta <= ta)
                        {
                            Eq.Add(qjID);
                            Ep.Add(pjID);
                            int[] e = new int[2];
                            e[0] = qjID;
                            e[1] = pjID;
                            E.Add(e);
                        }
                    }

                }

                if (n < E.Count)
                {
                    n = E.Count;
                    if (MatchedPair != null)
                    {
                        MatchedPair = E.ToList();
                        qRef = qiID;
                        pRef = piID;
                    }
                }

                Eq.Clear();
                Eq.TrimExcess();
                Eq = null;
                Ep.Clear();
                Ep.TrimExcess();
                Ep = null;
                E.Clear();
                E.TrimExcess();
                E = null;
            }

            return n;
        }

        public double sInv(mTriplet t, mTriplet r, ref mTriplet rMax)
        {
            double s = 1;

            double sDmax = Math.Abs(t.dMax - r.dMax);
            double sDmid = Math.Abs(t.dMid - r.dMid);
            double sDmin = Math.Abs(t.dMin - r.dMin);

            if (sDmax > tl || sDmid > tl || sDmin > tl)
            {
                s = 0;
            }
            else
            {
                mTriplet shift_r = shift(r);
                mTriplet shift2_r = shift(shift_r);
                double[] sParts = new double[3];
                sParts[0] = sPart(t, r);
                sParts[1] = sPart(t, shift_r);
                sParts[2] = sPart(t, shift2_r);

                s = sParts.Max();

                rMax = s == sParts[1] ? shift_r : s == sParts[2] ? shift2_r : r;
            }

            return s;
        }

        public mTriplet shift(mTriplet r)
        {
            mTriplet shifted = new mTriplet();

            shifted.pi[0] = r.pi[1];
            shifted.pi[1] = r.pi[2];
            shifted.pi[2] = r.pi[0];

            calIntrinsicFeatures(ref shifted);

            return shifted;
        }

        public void calIntrinsicFeatures(ref mTriplet r)
        {
            r.di[0] = EuclideanDistanceOf2Points(r.pi[0].p, r.pi[1].p);
            r.di[1] = EuclideanDistanceOf2Points(r.pi[0].p, r.pi[2].p);
            r.di[2] = EuclideanDistanceOf2Points(r.pi[1].p, r.pi[2].p);

            double dMax = r.di.Max();
            double dMin = r.di.Min();
            r.dMax = dMax;
            r.dMin = dMin;
            r.dMid = r.di.Where(d => d != dMax && d != dMin).FirstOrDefault();

            r.alphai[0] = ad2PI(r.pi[0].direction, ang(r.pi[1].p, r.pi[0].p));//
            r.alphai[1] = ad2PI(r.pi[0].direction, ang(r.pi[2].p, r.pi[0].p));
            r.alphai[2] = ad2PI(r.pi[1].direction, ang(r.pi[0].p, r.pi[1].p));
            r.alphai[3] = ad2PI(r.pi[1].direction, ang(r.pi[2].p, r.pi[1].p));//
            r.alphai[4] = ad2PI(r.pi[2].direction, ang(r.pi[0].p, r.pi[2].p));//
            r.alphai[5] = ad2PI(r.pi[2].direction, ang(r.pi[1].p, r.pi[2].p));

            r.betai[0] = ad2PI(r.pi[0].direction, r.pi[1].direction);
            r.betai[1] = ad2PI(r.pi[1].direction, r.pi[2].direction);
            r.betai[2] = ad2PI(r.pi[2].direction, r.pi[0].direction);
        }

        public double EuclideanDistanceOf2Points(Point a, Point b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public double EuclideanDistanceOf2Points(PointF a, PointF b)
        {
            return Math.Sqrt(Math.Pow((a.X - b.X), 2) + Math.Pow(a.Y - b.Y, 2));
        }

        /// <summary>
        /// adπ(alpha, beta) computes the minimum angle required to superpose two vectors with the same origin and angles alpha and beta respectively
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <returns>min(|alpha − beta|, 2π − |alpha − beta|)</returns>
        public double adPI(double alpha, double beta)
        {
            double absResult = Math.Abs(alpha - beta);
            double _2PI_AbsResult = (2 * Math.PI) - absResult;

            return absResult < _2PI_AbsResult ? absResult : _2PI_AbsResult;
        }

        /// <summary>
        /// ad2π(alpha, beta) computes the angle required to rotate a vector with angle beta in clockwise sense 
        /// to superpose it to another vector with the same origin and angle alpha
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        /// <returns>
        /// (beta - alpha) if beta > alpha
        /// (beta - alpha) + (2 * Math.PI) otherwise
        /// </returns>
        public double ad2PI(double alpha, double beta)
        {
            return beta > alpha ? (beta - alpha) : (beta - alpha) + (2 * Math.PI);
        }

        /// <summary>
        /// ang(pi, pj) computes the angle of the vector with initial point at pi and terminal point at pj
        /// </summary>
        /// <param name="pi"></param>
        /// <param name="pj"></param>
        /// <returns></returns>
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

        public double sPart(mTriplet t, mTriplet r)
        {
            double s = 1;

            double stheta = sTheta(t, r);
            double sdistance = sDistance(t, r);
            double salpha = sAlpha(t, r);
            double sbeta = sBeta(t, r);

            if (stheta == 0 || sdistance == 0 || salpha == 0 || sbeta == 0)
            {
                s = 0;
            }
            else
            {
                double sAll = (1 - sdistance) * (1 - salpha) * (1 - sbeta);
                s = 1 - sAll;
            }

            return s;
        }

        public double sTheta(mTriplet t, mTriplet r)
        {
            double s = 1;

            for (int i = 0; i < t.pi.Length; i++)
            {
                if (adPI(t.pi[i].direction, r.pi[i].direction) > (Math.PI / 4))
                {
                    s = 0;
                    break;
                }
            }

            return s;
        }

        public double sDistance(mTriplet t, mTriplet r)
        {
            double s = 1;
            double[] sdi = new double[t.di.Length];

            for (int i = 0; i < t.di.Length; i++)
            {
                sdi[i] = Math.Abs(t.di[i] - r.di[i]);
                if (sdi[i] > tl)
                {
                    s = 0;
                    break;
                }
            }

            if (s > 0)
            {
                s = 1 - (sdi.Max() / tl);
            }

            return s;
        }

        public double sAlpha(mTriplet t, mTriplet r)
        {
            double s = 1;
            double[] sAlphai = new double[t.alphai.Length];

            for (int i = 0; i < t.alphai.Length; i++)
            {
                sAlphai[i] = adPI(t.alphai[i], r.alphai[i]);
                if (sAlphai[i] > ta)
                {
                    s = 0;
                    break;
                }
            }

            if (s > 0)
            {
                s = 1 - (sAlphai.Max() / ta);
            }

            return s;
        }

        public double sBeta(mTriplet t, mTriplet r)
        {
            double s = 1;
            double[] sBetai = new double[t.betai.Length];

            for (int i = 0; i < t.betai.Length; i++)
            {
                sBetai[i] = adPI(t.betai[i], r.betai[i]);
                if (sBetai[i] > ta)
                {
                    s = 0;
                    break;
                }
            }

            if (s > 0)
            {
                s = 1 - (sBetai.Max() / ta);
            }

            return s;
        }
    }

    class mTripletPair
    {
        public mTriplet r;
        public mTriplet t;
        public double similarity;
    }
}
