# pragma once
#include "utils.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverTD{
public:
	struct Line
	{
		double k, b;
		Line():k(0),b(0){}
		void set(const double& _k, const double &_b){k=_k,b=_b;}
		Line(const double& k,const double& b): k(k), b(b) {}
		bool operator<(const Line &l) const {
			return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
		}
		bool operator==(const Line &l) const {return k == l.k;}
	};

	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1,
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						Array2d& colTime,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv = Array2d(0,DeltaT)) {
		auto ptPos1 = CpPos1.divideBezierPatch(divUvB1);
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		return primitiveCheck(ptPos1, ptVel1, ptPos2, ptVel2, colTime, bb, initTimeIntv);
	}

	static bool primitiveCheck(
			const std::array<Vector3d, ParamObj1::cntCp>& ptPos1,
			const std::array<Vector3d, ParamObj1::cntCp>& ptVel1,
			const std::array<Vector3d, ParamObj2::cntCp>& ptPos2,
			const std::array<Vector3d, ParamObj2::cntCp>& ptVel2,
			Array2d& colTime, const BoundingBoxType& bb,
			const Array2d& initTimeIntv,
			const double competingTime = INFT, bool* yielded = nullptr) {
		// yielded means that only partial axes were evaluated
		if(yielded != nullptr) *yielded = false;
		// Enlarge the time interval by a small margin so that the end points are correctly treated
		Array2d timeIntv(initTimeIntv[0]-1e-6, initTimeIntv[1]+1e-6);

		constexpr int MaxAxes = 15;
		std::array<Vector3d, MaxAxes> axes;
		const int axisCount = setAxes<ParamObj1, ParamObj2>(
			ptPos1, ptVel1, ptPos2, ptVel2, axes, bb, initTimeIntv[0]);

		std::array<Array2d, MaxAxes * 2> feasibleIntvs;
		int feasibleCount = 0;
		double minT = initTimeIntv[0], maxT = initTimeIntv[1];

		auto AxisCheck=[&](auto lines1, auto lines2){
			decltype(lines1) ch1;
			decltype(lines2) ch2;
			int chCount1 = 0, chCount2 = 0;
			for(auto& l:lines1)l.b += 1e-12 * std::abs(l.b);
			for(auto& l:lines2)l.b -= 1e-12 * std::abs(l.b);
			calcBoundaries(lines1, ch1, chCount1, true, timeIntv);
			calcBoundaries(lines2, ch2, chCount2, false, timeIntv);
			const auto intvT = boundaryIntersect(ch1, chCount1, ch2, chCount2, timeIntv);
			if(intvT[0]!=-1) feasibleIntvs[feasibleCount++] = intvT;
		};

		for(int axisId = 0; axisId < axisCount; ++axisId){
			const Vector3d& axis = axes[axisId];
			std::array<Line, ParamObj1::cntCp> ptLines1;
			std::array<Line, ParamObj2::cntCp> ptLines2;
			for(int i = 0; i < ParamObj1::cntCp; ++i)
				ptLines1[i] = Line(ptVel1[i].dot(axis), ptPos1[i].dot(axis));
			for(int i = 0; i < ParamObj2::cntCp; ++i)
				ptLines2[i] = Line(ptVel2[i].dot(axis), ptPos2[i].dot(axis));
			std::sort(ptLines1.begin(), ptLines1.end());
			std::sort(ptLines2.begin(), ptLines2.end());
			AxisCheck(ptLines1, ptLines2);
			AxisCheck(ptLines2, ptLines1);

			// Merge separating intervals connected to either end of the current candidate interval.
			double prevMinT, prevMaxT;
			do {
				prevMinT = minT;
				prevMaxT = maxT;
				for(int i = 0; i < feasibleCount; ++i){
					if(feasibleIntvs[i][0] <= minT)
						minT = std::max(minT, feasibleIntvs[i][1]);
					if(feasibleIntvs[i][1] >= maxT)
						maxT = std::min(maxT, feasibleIntvs[i][0]);
				}
			} while(minT != prevMinT || maxT != prevMaxT);

			// Separation has been certified over the entire candidate interval.
			if(minT > maxT) { colTime = Array2d(-1,-1); return false; }
			// minT can only increase as more axes are tested. Once it exceeds the best queued
			// lower bound, defer the remaining axes until this pair becomes competitive again.
			if(minT < maxT && minT > competingTime) {
				colTime = Array2d(minT, maxT);
				if(yielded != nullptr) *yielded = true;
				return true;
			}
		}

		colTime = Array2d(minT, maxT);
		return true;
	}

	template<std::size_t Count>
	static void calcBoundaries(std::array<Line, Count>& lines,
				std::array<Line, Count>& ch, int& chCount,
				const bool getMaxCH, const Array2d& tIntv) {
		// True for calculating the max func, false for calculating the min func
		if(!getMaxCH) std::reverse(lines.begin(), lines.end());
		const int lineCount = std::unique(lines.begin(), lines.end()) - lines.begin();

		chCount = 0;
		ch[chCount++] = lines[0];
		int alpha = 1;
		while(alpha < lineCount){
			int beta = chCount - 1;
			while(beta > 0){
				double chfp = (ch[beta].k-ch[beta-1].k)*(lines[alpha].b-ch[beta-1].b)
							-(lines[alpha].k-ch[beta-1].k)*(ch[beta].b-ch[beta-1].b);
				if(chfp>=0){
					--chCount;
					beta--;
				}
				else break;

			}
			if(beta==0){
				double chStart = tIntv[0]*(lines[alpha].k-ch[0].k)+(lines[alpha].b-ch[0].b);
				if((getMaxCH&&chStart>=0)||(!getMaxCH&&chStart<=0))
					--chCount;
			}
			if(chCount == 0)ch[chCount++] = lines[alpha];
			else{
				beta = chCount - 1;
				double chEnd = tIntv[1]*(lines[alpha].k-ch[beta].k)+(lines[alpha].b-ch[beta].b);
				if((getMaxCH&&chEnd>0)||(!getMaxCH&&chEnd<0))
					ch[chCount++] = lines[alpha];
			}
			alpha++;
		}
	}

	template<std::size_t Count1, std::size_t Count2>
	static Array2d boundaryIntersect(const std::array<Line, Count1>& ch1,
				const int chCount1, const std::array<Line, Count2>& ch2,
				const int chCount2, const Array2d& tIntv) {
		int id1=0, id2=0;

		// Find the left end of the intersection of the boundaries
		double intvL=-1, intvR=-1;
		if(ch1[0].k*tIntv[0]+ch1[0].b<ch2[0].k*tIntv[0]+ch2[0].b)intvL=tIntv[0];
		else{
			while(id1<chCount1&&id2<chCount2){
				if(ch1[id1].k>=ch2[id2].k){
					break;
				}
				double hifp1, hifp2;
				if(id1<chCount1-1)
					hifp1=(ch1[id1+1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1+1].b-ch2[id2].b);
				else
					hifp1=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2<chCount2-1)
					hifp2=(ch1[id1].k-ch2[id2+1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2+1].b);
				else
					hifp2=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(hifp1<0){
					if(hifp2<0){
						intvL = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						break;
					}
					else id2++;
				}
				else {
					id1++;
					if(hifp2<0);
					else id2++;
				}
			}
			if(intvL==-1||intvL>=tIntv[1])return Array2d(-1,-1);
		}

		// Find the right end of the intersection of the boundaries
		id1 = chCount1-1, id2 = chCount2-1;
		if((ch1[id1].k-ch2[id2].k)*tIntv[1]+(ch1[id1].b-ch2[id2].b)<0)intvR=tIntv[1];
		else{
			while(id1>=0&&id2>=0){
				if(ch1[id1].k<=ch2[id2].k){
					// std::cerr<<"WARNING: calculation of boundary intersection ends at strange slopes! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
					// exit(-1);
					return Array2d(-1,-1);
				}
				double hifp1, hifp2;
				if(id1>0)
					hifp1=(ch1[id1].k-ch2[id2].k)*(ch1[id1-1].b-ch2[id2].b)
							-(ch1[id1-1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b);
				else
					hifp1=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2>0)
					hifp2=(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2-1].b)
							-(ch1[id1].k-ch2[id2-1].k)*(ch1[id1].b-ch2[id2].b);
				else
					hifp2=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(hifp1<0){
					if(hifp2<0){
						intvR = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						break;
					}
					else id2--;
				}
				else {
					id1--;
					if(hifp2<0);
					else id2--;
				}
			}
			if(intvR==-1){
				std::cerr<<"WARNING: find the left end of boundary intersection but no right end! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
				std::cerr<<"("<<intvL<<", "<<intvR<<"), in candidate ("<<tIntv[0]<<", "<<tIntv[1]<<")\n";
				exit(-1);
			}
			if(intvR<=intvL)
				return Array2d(-1,-1);
		}

		if(intvL>intvR||intvL<tIntv[0]||intvR>tIntv[1]){
			// std::cerr<<"WARNING: you got an in-valid sub-time interval! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
			// std::cerr<<"("<<intvL<<", "<<intvR<<"), in candidate ("<<tIntv[0]<<", "<<tIntv[1]<<")\n";
			// exit(-1);
			return Array2d(-1,-1);
		}
		else return Array2d(intvL,intvR);
	}
public:
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1,
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2,
						const BoundingBoxType& bb,
						const double deltaDist,
						const double upperTime = DeltaT) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv;
			bool unresolved; // primitiveCheck yielded before testing all axes
			PatchPair(const ParamBound1& c1, const ParamBound2& c2,
					const Array2d& t = Array2d(0,DeltaT),
					const bool pending = false):
					pb1(c1), pb2(c2), tIntv(t), unresolved(pending) {}
			bool operator<(PatchPair const &o) const { return tIntv[0] > o.tIntv[0]; }
			double calcWidth() const{
				const double w1 = pb1.width(), w2 = pb2.width();
				return std::max(std::max(w1, w2), tIntv[1]-tIntv[0]);
			}
		};

		// Patch pairs compete by their current TOI lower bounds.
		std::priority_queue<PatchPair> heap;
		ParamBound1 initParam1;
		ParamBound2 initParam2;
		Array2d initTimeIntv(0,upperTime), colTime;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, bb, initTimeIntv))
			heap.emplace(initParam1, initParam2, colTime);
		while (!heap.empty()) {
			auto cur = heap.top();
			heap.pop();
			if(cur.unresolved){
				const double competingTime =
					bb == BoundingBoxType::OBB && !heap.empty() ? heap.top().tIntv[0] : INFT;
				bool yielded = false;
				if(!primitiveCheck(
						CpPos1.divideBezierPatch(cur.pb1),
						CpVel1.divideBezierPatch(cur.pb1),
						CpPos2.divideBezierPatch(cur.pb2),
						CpVel2.divideBezierPatch(cur.pb2),
						colTime, bb, cur.tIntv, competingTime, &yielded)) continue;
				cur.tIntv = colTime;
				cur.unresolved = yielded;
				if(yielded){
					heap.push(cur);
					continue;
				}
			}

			// Meets the precision requirement
			if (cur.calcWidth() < deltaDist) {
				uv1 = cur.pb1.centerParam();
				uv2 = cur.pb2.centerParam();
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			std::array<ParamBound1, 4> childBounds1;
			std::array<ParamBound2, 4> childBounds2;
			std::array<std::array<Vector3d, ParamObj1::cntCp>, 4> childPos1;
			std::array<std::array<Vector3d, ParamObj1::cntCp>, 4> childVel1;
			std::array<std::array<Vector3d, ParamObj2::cntCp>, 4> childPos2;
			std::array<std::array<Vector3d, ParamObj2::cntCp>, 4> childVel2;
			for(int i = 0; i < 4; ++i){
				childBounds1[i] = cur.pb1.interpSubpatchParam(i);
				childPos1[i] = CpPos1.divideBezierPatch(childBounds1[i]);
				childVel1[i] = CpVel1.divideBezierPatch(childBounds1[i]);
				childBounds2[i] = cur.pb2.interpSubpatchParam(i);
				childPos2[i] = CpPos2.divideBezierPatch(childBounds2[i]);
				childVel2[i] = CpVel2.divideBezierPatch(childBounds2[i]);
			}
			for (int i = 0; i < 4; i++) {
				const ParamBound1& divUvB1 = childBounds1[i];
				for (int j = 0; j < 4; j++) {
					const ParamBound2& divUvB2 = childBounds2[j];
					// Competitive yielding pays off for 15-axis OBB, but not for 3-axis AABB.
					const double competingTime =
						bb == BoundingBoxType::OBB && !heap.empty() ? heap.top().tIntv[0] : INFT;
					bool yielded = false;
					if (primitiveCheck(
							childPos1[i], childVel1[i], childPos2[j], childVel2[j],
							colTime, bb, cur.tIntv, competingTime, &yielded)){
						heap.emplace(divUvB1, divUvB2, colTime, yielded);
					}
				}
			}
		}

		return -1;
	}
};
