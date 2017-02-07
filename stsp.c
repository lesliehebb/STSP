#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define QUIET 0				//0 -> prints things, 1 -> only prints errors
#define QUIETMCMC 0			//0 -> mcmc prints as it goes, 1-> stays quiet (0 overridden by QUIET)
#define ALWAYSPRINTVIS 0
#define WHICHPRINTVIS 1		// what sort of visualization file to output (j-1,g-2)
#define ANYPRINTVIS 0		//0 for speed, overrides other PRINTVIS preferences
#define UNNORMALIZEOUTPUT 1	//1 -> undoes internal normalization before outputting 
#define PRINTEACHSPOT 0
#define DEFAULTFILENAME "window063_area.in"

#define ORDERSPOTSMETHOD 0			//0 ->order by latitude, 1->order by longitude, 2->hybrid ordering
#define DOWNFROMMAX 11				//how many light curve points down from the max to call the max (for normalization) 
#define MCMCTRACKMEM 1000000		//0 -> write mcmc tracker info to file, N -> buffer N outputs before writing
#define PRECALCPLANET 1				//1 -> optimize by precalculating planet effect at each time
#define COMBINEONESPOT	0			//whether to combine one spot at a time or a whole chain (should be 0)
#define COMBINERSINTHETA 0			//whether to use r or r*sin(theta) as the combination variable
#define USEOLDLDRING 1				//whether to use the old uniformly spaced limb darkening ring system

#define ALWAYSAVERAGETIME 0			//whether to average time when using maxsteps (in mcmc)
#define PRINTWHICHSPOT 1			//whether to print WHICHSPOT for final output (limits number of spots to 16 right now)
#define PRINTPLANETSPOTOVERLAP 0	//whether to print planet-spot overlap (lcgen only), may not work with flattened lightcurve

#define MAXSPOTS 10
#define MAXPLANETS 1
#define TWICEMAXSPOTSMAXPLANETS 20
#define FOURMAXSPOTSMAXPLANETS 40

#define MAXNLDRINGS 100

#define RND (((double)rand())/(double)RAND_MAX)
#define PIo2 1.5707963267948966192313216916398
#define PI 3.1415926535897932384626433832795
#define PIt2 6.283185307179586476925286766559
#define thPIo2 4.7123889803846898576939650749193
#define PIt2oday 0.000072722052166430399038487115353692
#define TOOBIG (1000000000)
#define TOOSMALL (-1000000000)
#define ACCEPTABLEERROR (0.000001)
#define ORDERSPOTNEARNESS (0.1745)

#define DEBUGMCMC 0
#define XYZDETAILS 1

int NLDRINGS;				//number of rings for limb darkening
char PRINTVIS;				//whether to print a vis file
char FIXTHETAS;				//whether to fix the latitudes of the spots
char USEDOWNFROMMAX;		//whether to use the DOWNFROMMAX for normalization (max brightness not given)
char CALCBRIGHTNESSFACTOR;	//whether to match normalization by calculating the brightness factor for each model
char FLATTENMODEL;			//whether to flatten the model outside of transits
char FIXSEEDEDONLYPHI;		//for partially seeded, fix only the phi's of seeds (not all parameters)
int  numseeded;
int numspots,numplanets;
long int starttime;
double spotfracdark;
double sigmaradius,sigmaangle;
FILE *outv,*outerr;
int errorflag;
int debuglcni,debugringi,debugtrial;
int pvasc,pvisc;
double *spotreport,*ldspotreport;
char whichspots[16];
unsigned int ttt[16]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};

#if DEBUGMCMC
FILE *dbmcmc;
#endif
#if XYZDETAILS
FILE *xyzdetail;
#endif

//FILE *debugellipse;

typedef struct 
{		//for ellipse-circle intersections
	double v,w;		//coords of intersection
	double gamma;	//angle from v-axis
}vwpoint;

typedef struct 
{
	int spot;		//spot it comes from
	char active;	//use to find blocked light if 1
	double centery,centerz;	//star frame coords
	double dstarcenter;	//distance from star center to ellipse center
	double smajor,sminor,area;	//shape of ellipse (s=semi)
	double fracdarkness;	//fractional darkness of ellipse (1=all dark)
	double vhaty,vhatz,whaty,whatz;	//ellipse frame unit vectors
	double mindcen,maxdcen;	//min and max distance from star center
	int numstarint;	//number of star-ellipse intersection points (-1 -> not yet calculated)
	vwpoint se[4];	//star-ellipse intersection points in elliipse frame
	int numcircleint[MAXPLANETS];	//number of intersections with each circle
	vwpoint circleint[MAXPLANETS][4];	//angle of planet intersection in ellipse frame, inside cricle from even to odd

}tellipse;

typedef struct 
{
	double centery,centerz;	//star frame coords
	double radius,area,rsq;	//shape of circle
	double mindcen,maxdcen;	//min and max distance from star center
	int numstarint;	//number of star-circle intersection points (-1 -> not yet calculated)
	double starintbeta[2];	//angle of star-circle intersections in star frame
	double starintdang;	//angle in planet frame (center) between star intersectoin points
	double scy[2],scz[2];	//coords of star-circle intersection in star frame
}tcircle;

typedef struct 
{
	char active;	//use to find blocked light if 1
	int ellipsen;	//ellipse it borders
	int ellipsearcn,scirclearcn;	//borders
	double area;	//total area
	double mindcen;	//min distance from star center
}tcresent;

typedef struct 
{
	int ellipsen;	//ellipse of which it is part
	vwpoint end[2];	//end points
}tellipsearc;

typedef struct 	//arc of real star's edge
{
	int ellipsen;	//ellipse it touches
	double radius;	//radius of arc
	vwpoint end[2];	//end points (ellipse frame)
	double ebeta[2];	//angles to end points in star frame
}tscirclearc;

typedef struct //arc of planet's edge
{
	int ellipsen;	//ellipse it intersects
	int circlen;	//circle it is part of
	vwpoint end[2];	//end points (ellipse frame)
	double edang;	//	angle between end points in circle frame
}tpcirclearc;							

typedef struct 
{
	int ellipsen;	//for coords
	vwpoint end[2];	//end points (ellipse frame) 0->1 ellipse on left
}tlinesegment;

typedef struct 	//sections are bulges bounded by a segment
{
	int ellipsearcn,linesegmentn;	//borders
	double area;	//total area
}tellipsesection;

typedef struct 
{
	int pcirclearcn,linesegmentn;	//borders
	double area;	//total area
}tcirclesection;

typedef struct 
{
	char cw;	//for sections with 1 ellipsearc and 1 scirclearc, cw=1 -> clockwise side of the cresent, cw=0 -> anticlockwise, cw=(-1) if undefined 
	int numellipsearcs,numscirclearcs;	//number of ellipse and scircle arcs in border
	int ellipsearcn[2],scirclearcn[2],linesegmentn;	//borders
	double area;	//total area
}tcresentsection;

typedef struct 	//ellipsesection - circlesection
{
	int ellipsen;	//number of ellipse
	int ellipsesectionn,circlesectionn;	//number of ellipse and circle section
	double area;	//total area
}tellipsepart;

typedef struct 	//cresentsection - circlesection
{
	int ellipsen;	//number of ellipse
	int cresentsectionn,circlesectionn;	//number of cresent and circle section
	double area;	//total area
}tcresentpart;

typedef struct 	//ellipse - circle
{
	int ellipsen,circlen;	//number of ellipse and circle
	double area;
}tellipsehole;

typedef struct 	//cresent - circle
{
	int cresentn,circlen;	//number of cresent and circle
	double area;
}tcresenthole;

typedef struct 	//number of each type that are being used
{
	int ellipse,circle,cresent,ellipsepart,cresentpart,ellipsehole,cresenthole;
	int ellipsesection,circlesection,cresentsection,ellipsearc,scirclearc,pcirclearc,linesegment;
}ttotalnum;

tellipse ellipse[MAXSPOTS];
tcircle circle[MAXPLANETS];
tcresent cresent[MAXSPOTS];
tellipsesection ellipsesection[TWICEMAXSPOTSMAXPLANETS];
tcirclesection circlesection[TWICEMAXSPOTSMAXPLANETS];
tcresentsection cresentsection[TWICEMAXSPOTSMAXPLANETS];
tellipsepart ellipsepart[TWICEMAXSPOTSMAXPLANETS];
tcresentpart cresentpart[TWICEMAXSPOTSMAXPLANETS];
tellipsehole ellipsehole[MAXSPOTS];
tcresenthole cresenthole[MAXSPOTS];
tellipsearc ellipsearc[FOURMAXSPOTSMAXPLANETS];
tscirclearc scirclearc[FOURMAXSPOTSMAXPLANETS];
tpcirclearc pcirclearc[TWICEMAXSPOTSMAXPLANETS];
tlinesegment linesegment[TWICEMAXSPOTSMAXPLANETS];
ttotalnum totalnum;
#if PRECALCPLANET
	tcircle *indcircle[MAXPLANETS];
	int *totalnumcircle;
#endif

typedef struct 
{
	double r;
	double fracdarkness;
	double alpha;  //central half angle of spot
	double flatdcen;	//distance to flat center of spot from star center
	double thetarel;  //theta relative to star rotation axis
	double poscoeff[5];  //for finding position at time t
	double theta;
	double phirel; //relative to star
	double phi;
	double psi;  //angle from x-axis
	double cospsi;
	double y,z;
	double area;  //unprojected area
}spotdata;

typedef struct 
{
	double r;
	double rsq;  //radius squared
	double area;
	double omega;
	double theta;  //tilt of rotation axis towards us
	double phi;  //rotation at particular time
	double dintensity[MAXNLDRINGS];  //ring k is dintensity[k] brighter than ring k+1
	double ringr[MAXNLDRINGS];  //outer radius of ring k
	double maxlight;	//from adding up limb darkening rings
	double brightnessfactor;	//to get brightness of unoccluded star (brighter than max of data)
}stardata;

typedef struct 
{
	double r;
	double rsq;					//radius squared
	double area;
	double omegaorbit;			//angular speed of orbit rad/s
	double orbitsemimajor;		//semimajor axis length
	double eccentricity;	
	double thetaorbit,phiorbit;		//theta and phi of the orbital axis
	double orbitangleomega;			//angle twixt line of nodes and semimajor axis (towards perastrion)
	double tmiddleoftransit;		//t in planetpos fuction at which planet is in middle of eclipse, to match with t=0 of data
	double xhat[3],yhat[3];			//unit vectors for orbital ellipse frame
	double lct0;					//at t=0, light curve time=lct0 (used for data) 
	double x,y,z;						//coords at particular time
}planetdata;


long int timecheck(void)
{
	time_t time0;

	time0=time(NULL);
	return (long int)time0;
}
char isanglebetween(double cw,double ccw,double between)
{
	double diftot,difq;
	if(cw>TOOBIG||ccw>TOOBIG||cw<TOOSMALL||ccw<TOOSMALL)
	{
		errorflag=1;
		return 0;
	}
	while(cw<0) cw+=PIt2;
	while(ccw<0) ccw+=PIt2;
	while(between<0) between+=PIt2;
	diftot=ccw-cw;
	if(diftot<0) diftot+=PIt2;
	difq=between-cw;
	if(difq<0) difq+=PIt2;
	if(difq<diftot)
		return 1;
	else
		return 0;
}
void vwtoyz(double v,double w,int en,double *y,double *z)
{
	(*y)=ellipse[en].centery+v*ellipse[en].vhaty+w*ellipse[en].whaty;
	(*z)=ellipse[en].centerz+v*ellipse[en].vhatz+w*ellipse[en].whatz;
}
double vwpointtobeta(vwpoint p,int en)
{
	double y,z,beta;
	vwtoyz(p.v,p.w,en,&y,&z);
	beta=atan2(z,y);
	if(beta<0)	beta+=PIt2;
	return beta;
}
double sebeta(int sen,int en)
{
	double y,z,beta;
	vwtoyz(ellipse[en].se[sen].v,ellipse[en].se[sen].w,en,&y,&z);
	beta=atan2(z,y);
	if(beta<0)	beta+=PIt2;
	return beta;
}
char wonleft(int ln,double w0)
{
	if((linesegment[ln].end[1].v-linesegment[ln].end[0].v)*(w0-linesegment[ln].end[0].w)+linesegment[ln].end[0].v*(linesegment[ln].end[1].w-linesegment[ln].end[0].w)>0)
		return 1;
	else
		return 0;
}
int ellipsecircleintersection(double a,double b,double v0,double w0,double rsq,vwpoint p[4])
{	//intersection of ellipse (major axis 'a' along v-axis, minor axis 'b' along w-axis) and circle at v0,w0 of radius squared rsq
	int i;				//ellipse inside circle from even p to odd p
	double g,c00,c0,c1,vdif,wdif,gpos,gneg,gmid;
	double gstep=0.01,gmax;
	gmax=PIt2-gstep;
	
	for(i=0;i<4;i++)
		p[i].gamma=(-1);  

	vdif=a-v0;
	c00=vdif*vdif+w0*w0-rsq; //g of 0.00
	c1=c00;
	for(g=gstep;g<=PIt2;g+=gstep)
	{
		c0=c1;
		vdif=a*cos(g)-v0;
		wdif=b*sin(g)-w0;
		c1=vdif*vdif+wdif*wdif-rsq;  //squared distance of point on ellipse from center of circle minus radius squared
		if(g>gmax) //very careful to go all the way round
			c1=c00;
		if((c0>=0.0&&c1<0.0)||(c0<0.0&&c1>=0.0))  //change of sign
		{
			if(g<=gmax)
			{
				gpos=gneg=g;
				if(c0>=0.0)
					gpos-=gstep;
				else
					gneg-=gstep;
			}
			else
			{	//being carefull to match up at the other side
				if(c0>=0)
				{
					gpos=g-gstep;
					gneg=PIt2;
				}
				else
				{
					gpos=PIt2;
					gneg=g-gstep;
				}
				g=6.3; //to exit loop
			}
			for(i=0;i<8;i++)
			{
				gmid=(gneg+gpos)/2.0;
				vdif=a*cos(gmid)-v0;
				wdif=b*sin(gmid)-w0;
				if(vdif*vdif+wdif*wdif-rsq>=0)
					gpos=gmid;
				else
					gneg=gmid;
			}
			gmid=(gneg+gpos)/2.0;
			if(gpos<gneg)
				if(p[0].gamma==(-1))
					p[0].gamma=gmid;      //this 'gamma' is a parameter of the curve, not the 
				else if(p[2].gamma==(-1))		//actual angle, it's adjusted below, to be the angle
					p[2].gamma=gmid;
				else
					errorflag=2;
			else
				if(p[0].gamma==(-1))
					if(p[3].gamma==(-1))
						p[3].gamma=gmid;
					else
						errorflag=2;
				else
					if(p[1].gamma==(-1))
						p[1].gamma=gmid;
					else if(p[3].gamma==(-1))
						p[3].gamma=gmid;
					else
						errorflag=2;
		}
	}
	if(p[0].gamma==(-1))  //no points
		return 0;
	if(p[2].gamma==(-1))  //2 points
	{
		if(p[1].gamma==(-1))
			if(p[3].gamma==(-1))
				errorflag=2;
			else
				p[1].gamma=p[3].gamma;
		for(i=0;i<2;i++)
		{
			p[i].v=a*cos(p[i].gamma);
			p[i].w=b*sin(p[i].gamma);
			p[i].gamma=atan2(p[i].w,p[i].v);
			if(p[i].gamma<0)
				p[i].gamma+=PIt2;
		}
		return 2;
	}
	else	//4 points
	{
		for(i=0;i<4;i++)
		{
			if(p[i].gamma==(-1))
				errorflag=2;
			p[i].v=a*cos(p[i].gamma);
			p[i].w=b*sin(p[i].gamma);
			p[i].gamma=atan2(p[i].w,p[i].v);
			if(p[i].gamma<0)
				p[i].gamma+=PIt2;
		}
		return 4;
	}
}
int circlestarintersection(double r,double dcen,double rcir,double v[2],double w[2])
{
	//finds intersection points of star (radius r) with circle dcen from center, of radius rcir
	//only returns positive-v point of each pair, returns number of points
	w[0]=(r*r-dcen*dcen-rcir*rcir)/(2.0*dcen);
	if(w[0]>=rcir||w[0]<=(-rcir))
		return 0;
	v[0]=sqrt(rcir*rcir-w[0]*w[0]);
	return 1;
}
int ellipsestarintersection(double r,double dcen,double a,double b,double v[2],double w[2])
{	
	//finds intersection points of star (radius r) with ellipse dcen from center, semimajor axis a, semiminor axis b
	//only sets positive-v point of each pair, returns number of pairs of points (0,1,or 2)
	int k;
	double wp,wn,t0,t1,t2,aobsq;
	if(a==b)
		return circlestarintersection(r,dcen,a,v,w);
	aobsq=a*a/(b*b);
	t1=4.0*dcen*dcen-4.0*(1.0-aobsq)*(a*a+dcen*dcen-r*r);
	if(t1<0)
		return 0;
	t1=sqrt(t1);
	t0=(-2.0)*dcen;
	t2=2.0*(1.0-aobsq);
	wp=(t0+t1)/t2;
	wn=(t0-t1)/t2;
	k=0;
	if(wn>(-b)&&wn<b)		//if 2 pairs, v[0],w[0] will be the lower pair (negative w)
	{
		w[k]=wn;
		v[k]=sqrt(a*a-aobsq*wn*wn);
		k++;
	}
	if(wp>(-b)&&wp<b)
	{
		w[k]=wp;
		v[k]=sqrt(a*a-aobsq*wp*wp);
		k++;
	}
	if(k==2&&w[0]>w[1])
	{	//make sure
		t0=v[0]; t1=w[0];
		v[0]=v[1]; w[0]=w[1];
		v[1]=t0; w[1]=t1;
	}
	return k;
}
double areacirclesection(double r,double dang)
{	//area of a circular section radius r, from angle ang0 to ang1 (tangent (v hat) is ang=0) -bordered by arc and chord (not a sector)
	double alpha,halfchord;
	alpha=dang/2.0;
	if(alpha>PIo2)
	{
		alpha=PI-alpha;
		halfchord=r*sin(alpha);
		return (PI-alpha)*r*r+sqrt(r*r-halfchord*halfchord)*halfchord;
	}
	else
	{
		halfchord=r*sin(alpha);
		return alpha*r*r-sqrt(r*r-halfchord*halfchord)*halfchord;
	}
}
double areacirclesection2(double r,double ang0,double ang1)
{	//area of a circular section radius r, from angle ang0 to ang1 (tangent (v hat) is ang=0)
	double alpha,halfchord;
	if(ang0>ang1)
		ang1+=PIt2;
	alpha=(ang1-ang0)/2.0;
	if(alpha>PIo2)
	{
		alpha=PI-alpha;
		halfchord=r*sin(alpha);
		return (PI-alpha)*r*r+sqrt(r*r-halfchord*halfchord)*halfchord;
	}
	else
	{
		halfchord=r*sin(alpha);
		return alpha*r*r-sqrt(r*r-halfchord*halfchord)*halfchord;
	}
}
double areaellipsesection(double a,double b,double ang0,double ang1)
{	//	area of an elliptical section, major axis a (along v hat) minor axis b (along w hat),
	//	from angle ang0 to ang1 (tangent (v-hat) is ang=0)
	double f,v,w,cang0,cang1,ep;
	if(b<0)
		b=(-b);  //for ellipses round back
	f=a/b;  //stretch factor along w-hat
	ep=atan(a*tan(ang0)/b);  //from -PIo2 to PIo2
	if(ang0>PIo2&&ang0<thPIo2)
		ep+=PI;
	v=a*cos(ep);	
	w=b*sin(ep);
	w*=f;				//stretch to circle
	cang0=atan2(w,v);
	if(cang0<0)
		cang0+=PIt2;
	ep=atan(a*tan(ang1)/b);
	if(ang1>PIo2&&ang1<thPIo2)
		ep+=PI;
	v=a*cos(ep);
	w=b*sin(ep);
	w*=f;
	cang1=atan2(w,v);
	if(cang1<0)
		cang1+=PIt2;	//find area and compress back to ellipse
	return areacirclesection2(a,cang0,cang1)/f;	
}
double areatriangle(vwpoint a,vwpoint b,vwpoint c)
{			//area of triangle
	double qv[2],qw[2],area;

	qv[0]=b.v-a.v;
	qw[0]=b.w-a.w;
	qv[1]=c.v-a.v;
	qw[1]=c.w-a.w;
	area=qv[0]*qw[1]-qw[0]*qv[1];  //really just cross products
	return area/2.0;
}
double areaquad(vwpoint p[4])
{			//darkness of quadrilateral p0p1p2p3 given in v-w coords
	int i;
	double qv[3],qw[3],a0,a1;
	for(i=0;i<3;i++)
	{
		qv[i]=p[i].v-p[3].v;
		qw[i]=p[i].w-p[3].w;
	}
	a0=qv[0]*qw[1]-qw[0]*qv[1];  //really just cross products
	a1=qv[1]*qw[2]-qw[1]*qv[2];
	if(a0<0||a1<0)
		errorflag=11;  //negative area error
	return (a0+a1)/2.0;
}
void updatespots(spotdata spot[MAXSPOTS],double starradius,double starphi)
{
	int i;
	double cp,sp,tmp;
	for(i=0;i<numspots;i++)
	{
		sp=spot[i].phirel+starphi;
		cp=cos(sp);
		sp=sin(sp);
		tmp=spot[i].poscoeff[0]*cp+spot[i].poscoeff[1];
		if(tmp>(-1.0)&&tmp<1.0)
			spot[i].theta=acos(tmp); //should be [0-pi]
		else if(tmp>=1)
			spot[i].theta=0.0;
		else
			spot[i].theta=PI;
		tmp=(spot[i].poscoeff[2]*cp+spot[i].poscoeff[3])/sin(spot[i].theta);
		if(tmp>(-1.0)&&tmp<1.0)
			spot[i].phi=acos(tmp);
		else if(tmp>=1)
			spot[i].phi=0.0;
		else
			spot[i].phi=PI;
		if(spot[i].poscoeff[4]*sp/sin(spot[i].theta)<0.0)
			spot[i].phi=PIt2-spot[i].phi;
		spot[i].cospsi=sin(spot[i].theta)*cos(spot[i].phi);
		spot[i].psi=acos(spot[i].cospsi);
		spot[i].y=starradius*cos(spot[i].alpha)*sin(spot[i].theta)*sin(spot[i].phi);
		spot[i].z=starradius*cos(spot[i].alpha)*cos(spot[i].theta);
	}
}
void zerototalnums(void)
{
	totalnum.circle=0; totalnum.circlesection=0; totalnum.cresent=0; totalnum.cresenthole=0;
	totalnum.cresentpart=0; totalnum.cresentsection=0; totalnum.ellipse=0; totalnum.ellipsearc=0;
	totalnum.ellipsehole=0; totalnum.ellipsepart=0; totalnum.ellipsesection=0; totalnum.linesegment=0;
	totalnum.pcirclearc=0; totalnum.scirclearc=0;
}
void createcircle(planetdata oneplanet,double starradius,double starradiussq)
{
	double d,a,c,anga,angb;
	circle[totalnum.circle].area=oneplanet.area;
	circle[totalnum.circle].centery=oneplanet.y;
	circle[totalnum.circle].centerz=oneplanet.z;
	d=sqrt(oneplanet.y*oneplanet.y+oneplanet.z*oneplanet.z);
	circle[totalnum.circle].maxdcen=d+oneplanet.r;
	circle[totalnum.circle].mindcen=d-oneplanet.r;
	circle[totalnum.circle].radius=oneplanet.r;
	circle[totalnum.circle].rsq=oneplanet.rsq;
	if((d>starradius+oneplanet.r)||(d<oneplanet.r-starradius))	//check intersections with real star
		circle[totalnum.circle].numstarint=0;
	else
	{
		a=(starradiussq-oneplanet.rsq+d*d)/(2*d);
		c=sqrt(starradiussq-a*a);
		angb=atan2(oneplanet.z,oneplanet.y);
		anga=acos(a/starradius);
		circle[totalnum.circle].starintbeta[0]=angb-anga;
		if(circle[totalnum.circle].starintbeta[0]<0)	circle[totalnum.circle].starintbeta[0]+=PIt2;
		circle[totalnum.circle].starintbeta[1]=angb+anga;
		if(circle[totalnum.circle].starintbeta[1]<0)	circle[totalnum.circle].starintbeta[1]+=PIt2;
		a=asin(c/oneplanet.r);
		if(d*d<starradiussq-oneplanet.rsq)
			a=PI-a;	//mostly inside real star
		circle[totalnum.circle].starintdang=2.0*a;
		circle[totalnum.circle].numstarint=2;
	}
	totalnum.circle++;
}
void createellipsecresent(int spotn,spotdata onespot,double starradius)
{
	double d,x,y,z;
	ellipse[totalnum.ellipse].spot=spotn;
	ellipse[totalnum.ellipse].area=onespot.area*onespot.cospsi;
	ellipse[totalnum.ellipse].centery=onespot.y;
	ellipse[totalnum.ellipse].centerz=onespot.z;
	ellipse[totalnum.ellipse].fracdarkness=onespot.fracdarkness;
	ellipse[totalnum.ellipse].smajor=onespot.r;
	ellipse[totalnum.ellipse].sminor=onespot.r*onespot.cospsi;
	if(ellipse[totalnum.ellipse].sminor<0) ellipse[totalnum.ellipse].sminor=(-ellipse[totalnum.ellipse].sminor);
	d=sqrt(onespot.y*onespot.y+onespot.z*onespot.z);
	ellipse[totalnum.ellipse].maxdcen=sqrt(ellipse[totalnum.ellipse].smajor*ellipse[totalnum.ellipse].smajor+d*d*(1.0+1.0/(ellipse[totalnum.ellipse].smajor*ellipse[totalnum.ellipse].smajor/(ellipse[totalnum.ellipse].sminor*ellipse[totalnum.ellipse].sminor)-1.0)));
	ellipse[totalnum.ellipse].mindcen=d-ellipse[totalnum.ellipse].sminor;
	ellipse[totalnum.ellipse].dstarcenter=d;
	if(d>0)
	{
		ellipse[totalnum.ellipse].whaty=ellipse[totalnum.ellipse].centery/d;
		ellipse[totalnum.ellipse].whatz=ellipse[totalnum.ellipse].centerz/d;
		ellipse[totalnum.ellipse].vhaty=ellipse[totalnum.ellipse].whatz;
		ellipse[totalnum.ellipse].vhatz=(-ellipse[totalnum.ellipse].whaty);
	}
	else
	{
		ellipse[totalnum.ellipse].whaty=0.0; ellipse[totalnum.ellipse].whatz=1.0;
		ellipse[totalnum.ellipse].vhaty=1.0; ellipse[totalnum.ellipse].vhatz=0.0;
	}
	if(onespot.psi<PIo2)
		ellipse[totalnum.ellipse].active=1;
	else
		ellipse[totalnum.ellipse].active=0;
	if(onespot.psi+onespot.alpha<=PIo2)
		ellipse[totalnum.ellipse].numstarint=0; //for real star
	else	//create cresent
	{
		ellipse[totalnum.ellipse].numstarint=0;	//for real star, it's tangent -> no 'intersection'
		cresent[totalnum.cresent].mindcen=d+ellipse[totalnum.ellipse].sminor;
		cresent[totalnum.cresent].ellipsen=totalnum.ellipse;
		ellipsearc[totalnum.ellipsearc].ellipsen=totalnum.ellipse;
		if(onespot.psi<PIo2)
			x=onespot.flatdcen*tan(PIo2-onespot.psi);
		else
			x=onespot.flatdcen*tan(onespot.psi-PIo2);
		ellipsearc[totalnum.ellipsearc].end[0].v=sqrt(ellipse[totalnum.ellipse].smajor*ellipse[totalnum.ellipse].smajor-x*x);
		ellipsearc[totalnum.ellipsearc].end[0].w=x*onespot.cospsi;
		if(ellipsearc[totalnum.ellipsearc].end[0].w<0)	ellipsearc[totalnum.ellipsearc].end[0].w=(-ellipsearc[totalnum.ellipsearc].end[0].w);
		ellipsearc[totalnum.ellipsearc].end[0].gamma=atan2(ellipsearc[totalnum.ellipsearc].end[0].w,ellipsearc[totalnum.ellipsearc].end[0].v);
		ellipsearc[totalnum.ellipsearc].end[1].v=(-ellipsearc[totalnum.ellipsearc].end[0].v);
		ellipsearc[totalnum.ellipsearc].end[1].w=ellipsearc[totalnum.ellipsearc].end[0].w;
		ellipsearc[totalnum.ellipsearc].end[1].gamma=PI-ellipsearc[totalnum.ellipsearc].end[0].gamma;
		cresent[totalnum.cresent].ellipsearcn=totalnum.ellipsearc;
		scirclearc[totalnum.scirclearc].ellipsen=totalnum.ellipse;
		scirclearc[totalnum.scirclearc].end[0].v=ellipsearc[totalnum.ellipsearc].end[0].v;
		scirclearc[totalnum.scirclearc].end[0].w=ellipsearc[totalnum.ellipsearc].end[0].w;
		scirclearc[totalnum.scirclearc].end[0].gamma=ellipsearc[totalnum.ellipsearc].end[0].gamma;
		scirclearc[totalnum.scirclearc].end[1].v=ellipsearc[totalnum.ellipsearc].end[1].v;
		scirclearc[totalnum.scirclearc].end[1].w=ellipsearc[totalnum.ellipsearc].end[1].w;
		scirclearc[totalnum.scirclearc].end[1].gamma=ellipsearc[totalnum.ellipsearc].end[1].gamma;
		vwtoyz(scirclearc[totalnum.scirclearc].end[0].v,scirclearc[totalnum.scirclearc].end[0].w,totalnum.ellipse,&y,&z);
		scirclearc[totalnum.scirclearc].ebeta[0]=atan2(z,y);
		if(scirclearc[totalnum.scirclearc].ebeta[0]<0)	scirclearc[totalnum.scirclearc].ebeta[0]+=PIt2;
		vwtoyz(scirclearc[totalnum.scirclearc].end[1].v,scirclearc[totalnum.scirclearc].end[1].w,totalnum.ellipse,&y,&z);
		scirclearc[totalnum.scirclearc].ebeta[1]=atan2(z,y);
		if(scirclearc[totalnum.scirclearc].ebeta[1]<0)	scirclearc[totalnum.scirclearc].ebeta[1]+=PIt2;
		scirclearc[totalnum.scirclearc].radius=starradius;
		cresent[totalnum.cresent].scirclearcn=totalnum.scirclearc;
		cresent[totalnum.cresent].area=areacirclesection2(scirclearc[totalnum.scirclearc].radius,scirclearc[totalnum.scirclearc].ebeta[0],scirclearc[totalnum.scirclearc].ebeta[1])
			-areaellipsesection(ellipse[totalnum.ellipse].smajor,ellipse[totalnum.ellipse].sminor,ellipsearc[totalnum.ellipsearc].end[0].gamma,ellipsearc[totalnum.ellipsearc].end[1].gamma);
		cresent[totalnum.cresent].active=1;
		totalnum.ellipsearc++;
		totalnum.scirclearc++;
		totalnum.cresent++;
	}

	totalnum.ellipse++;
}
void createellipsepart(vwpoint ecw,vwpoint eccw,int en,int cn,double vc,double wc)
{
	double ccwa,cwa;
	vwpoint t;
	linesegment[totalnum.linesegment].ellipsen=en;
	linesegment[totalnum.linesegment].end[0].gamma=eccw.gamma;
	linesegment[totalnum.linesegment].end[0].v=eccw.v;
	linesegment[totalnum.linesegment].end[0].w=eccw.w;
	linesegment[totalnum.linesegment].end[1].gamma=ecw.gamma;
	linesegment[totalnum.linesegment].end[1].v=ecw.v;
	linesegment[totalnum.linesegment].end[1].w=ecw.w;

	ellipsearc[totalnum.ellipsearc].ellipsen=en;
	ellipsearc[totalnum.ellipsearc].end[0].gamma=ecw.gamma;
	ellipsearc[totalnum.ellipsearc].end[0].v=ecw.v;
	ellipsearc[totalnum.ellipsearc].end[0].w=ecw.w;
	ellipsearc[totalnum.ellipsearc].end[1].gamma=eccw.gamma;
	ellipsearc[totalnum.ellipsearc].end[1].v=eccw.v;
	ellipsearc[totalnum.ellipsearc].end[1].w=eccw.w;

	pcirclearc[totalnum.pcirclearc].circlen=cn;
	pcirclearc[totalnum.pcirclearc].ellipsen=en;
	pcirclearc[totalnum.pcirclearc].end[0].gamma=ecw.gamma;
	pcirclearc[totalnum.pcirclearc].end[0].v=ecw.v;
	pcirclearc[totalnum.pcirclearc].end[0].w=ecw.w;
	pcirclearc[totalnum.pcirclearc].end[1].gamma=eccw.gamma;
	pcirclearc[totalnum.pcirclearc].end[1].v=eccw.v;
	pcirclearc[totalnum.pcirclearc].end[1].w=eccw.w;
	cwa=atan2(ecw.w-wc,ecw.v-vc);
	if(cwa<0)	cwa+=PIt2;
	ccwa=atan2(eccw.w-wc,eccw.v-vc);
	if(ccwa<0)	ccwa+=PIt2;
	if(ccwa<cwa)
		ccwa+=PIt2;
	pcirclearc[totalnum.pcirclearc].edang=ccwa-cwa;

	ellipsesection[totalnum.ellipsesection].ellipsearcn=totalnum.ellipsearc;
	ellipsesection[totalnum.ellipsesection].linesegmentn=totalnum.linesegment;
	ellipsesection[totalnum.ellipsesection].area=areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ecw.gamma,eccw.gamma);

	circlesection[totalnum.circlesection].pcirclearcn=totalnum.pcirclearc;
	circlesection[totalnum.circlesection].linesegmentn=totalnum.linesegment;
	circlesection[totalnum.circlesection].area=areacirclesection(circle[cn].radius,pcirclearc[totalnum.pcirclearc].edang);

	ellipsepart[totalnum.ellipsepart].ellipsen=en;
	ellipsepart[totalnum.ellipsepart].ellipsesectionn=totalnum.ellipsesection;
	ellipsepart[totalnum.ellipsepart].circlesectionn=totalnum.circlesection;
	ellipsepart[totalnum.ellipsepart].area=ellipsesection[totalnum.ellipsesection].area-circlesection[totalnum.circlesection].area;
	if(ellipsepart[totalnum.ellipsepart].area<0)
		if(ellipse[en].area<ACCEPTABLEERROR)
		{
			t.v=pcirclearc[totalnum.pcirclearc].end[0].v; t.w=pcirclearc[totalnum.pcirclearc].end[0].w; t.gamma=pcirclearc[totalnum.pcirclearc].end[0].gamma;
			pcirclearc[totalnum.pcirclearc].end[0].v=pcirclearc[totalnum.pcirclearc].end[1].v; pcirclearc[totalnum.pcirclearc].end[0].w=pcirclearc[totalnum.pcirclearc].end[1].w; pcirclearc[totalnum.pcirclearc].end[0].gamma=pcirclearc[totalnum.pcirclearc].end[1].gamma;
			pcirclearc[totalnum.pcirclearc].end[1].v=t.v; pcirclearc[totalnum.pcirclearc].end[1].w=t.w; pcirclearc[totalnum.pcirclearc].end[1].gamma=t.gamma;
			pcirclearc[totalnum.pcirclearc].edang=PIt2-pcirclearc[totalnum.pcirclearc].edang;
			circlesection[totalnum.circlesection].area=areacirclesection(circle[cn].radius,pcirclearc[totalnum.pcirclearc].edang);
			ellipsepart[totalnum.ellipsepart].area=ellipsesection[totalnum.ellipsesection].area-circlesection[totalnum.circlesection].area;
		}
	if(ellipsepart[totalnum.ellipsepart].area<0)
		errorflag=4;
	totalnum.linesegment++;
	totalnum.ellipsearc++;
	totalnum.pcirclearc++;
	totalnum.ellipsesection++;
	totalnum.circlesection++;
	totalnum.ellipsepart++;
}
void createcresentpart20(vwpoint p0,vwpoint p1,int cren,int en,int cirn,double vc,double wc)
{
	double ccwa,cwa;

	linesegment[totalnum.linesegment].ellipsen=en;
	linesegment[totalnum.linesegment].end[0].v=p1.v;
	linesegment[totalnum.linesegment].end[0].w=p1.w;
	linesegment[totalnum.linesegment].end[0].gamma=p1.gamma;
	linesegment[totalnum.linesegment].end[1].v=p0.v;
	linesegment[totalnum.linesegment].end[1].w=p0.w;
	linesegment[totalnum.linesegment].end[1].gamma=p0.gamma;

	cresentsection[totalnum.cresentsection].numellipsearcs=2;
	cresentsection[totalnum.cresentsection].numscirclearcs=1;
	cresentsection[totalnum.cresentsection].linesegmentn=totalnum.linesegment;
	cresentsection[totalnum.cresentsection].cw=(-1);
	ellipsearc[totalnum.ellipsearc].ellipsen=en;
	ellipsearc[totalnum.ellipsearc].end[0].v=ellipsearc[cresent[cren].ellipsearcn].end[0].v;
	ellipsearc[totalnum.ellipsearc].end[0].w=ellipsearc[cresent[cren].ellipsearcn].end[0].w;
	ellipsearc[totalnum.ellipsearc].end[0].gamma=ellipsearc[cresent[cren].ellipsearcn].end[0].gamma;
	ellipsearc[totalnum.ellipsearc].end[1].v=p0.v;
	ellipsearc[totalnum.ellipsearc].end[1].w=p0.w;
	ellipsearc[totalnum.ellipsearc].end[1].gamma=p0.gamma;
	cresentsection[totalnum.cresentsection].ellipsearcn[0]=totalnum.ellipsearc;	//from cw end of ellipse arc to L0
	totalnum.ellipsearc++;
	ellipsearc[totalnum.ellipsearc].ellipsen=en;
	ellipsearc[totalnum.ellipsearc].end[1].v=ellipsearc[cresent[cren].ellipsearcn].end[1].v;
	ellipsearc[totalnum.ellipsearc].end[1].w=ellipsearc[cresent[cren].ellipsearcn].end[1].w;
	ellipsearc[totalnum.ellipsearc].end[1].gamma=ellipsearc[cresent[cren].ellipsearcn].end[1].gamma;
	ellipsearc[totalnum.ellipsearc].end[0].v=p1.v;
	ellipsearc[totalnum.ellipsearc].end[0].w=p1.w;
	ellipsearc[totalnum.ellipsearc].end[0].gamma=p1.gamma;
	cresentsection[totalnum.cresentsection].ellipsearcn[1]=totalnum.ellipsearc;	//from L1 to ccw end of ellipse arc
	totalnum.ellipsearc++;
	cresentsection[totalnum.cresentsection].scirclearcn[0]=cresent[cren].scirclearcn;
	cresentsection[totalnum.cresentsection].area=cresent[cren].area+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,p0.gamma,p1.gamma);

	circlesection[totalnum.circlesection].linesegmentn=totalnum.linesegment;
	pcirclearc[totalnum.pcirclearc].circlen=cirn;
	pcirclearc[totalnum.pcirclearc].ellipsen=en;
	pcirclearc[totalnum.pcirclearc].end[0].v=p0.v;
	pcirclearc[totalnum.pcirclearc].end[0].w=p0.w;
	pcirclearc[totalnum.pcirclearc].end[0].gamma=p0.gamma;
	pcirclearc[totalnum.pcirclearc].end[1].v=p1.v;
	pcirclearc[totalnum.pcirclearc].end[1].w=p1.w;
	pcirclearc[totalnum.pcirclearc].end[1].gamma=p1.gamma;
	cwa=atan2(pcirclearc[totalnum.pcirclearc].end[0].w-wc,pcirclearc[totalnum.pcirclearc].end[0].v-vc);
	if(cwa<0)	cwa+=PIt2;
	ccwa=atan2(pcirclearc[totalnum.pcirclearc].end[1].w-wc,pcirclearc[totalnum.pcirclearc].end[1].v-vc);
	if(ccwa<0)	ccwa+=PIt2;
	if(ccwa<cwa)
		ccwa+=PIt2;
	pcirclearc[totalnum.pcirclearc].edang=ccwa-cwa;
	circlesection[totalnum.circlesection].pcirclearcn=totalnum.pcirclearc;
	totalnum.pcirclearc++;
	circlesection[totalnum.circlesection].area=areacirclesection(circle[cirn].radius,pcirclearc[circlesection[totalnum.circlesection].pcirclearcn].edang);

	totalnum.linesegment++;
	cresentpart[totalnum.cresentpart].circlesectionn=totalnum.circlesection;
	totalnum.circlesection++;
	cresentpart[totalnum.cresentpart].cresentsectionn=totalnum.cresentsection;
	totalnum.cresentsection++;
	cresentpart[totalnum.cresentpart].area=cresentsection[cresentpart[totalnum.cresentpart].cresentsectionn].area-circlesection[cresentpart[totalnum.cresentpart].circlesectionn].area;
	cresentpart[totalnum.cresentpart].ellipsen=en;
	totalnum.cresentpart++;
}
void createcresentpart02(vwpoint q0,vwpoint q1,int cren,int en,int cirn,double vc,double wc,double starradius)
{
	double ccwa,cwa;

	linesegment[totalnum.linesegment].ellipsen=en;
	linesegment[totalnum.linesegment].end[0].v=q0.v;
	linesegment[totalnum.linesegment].end[0].w=q0.w;
	linesegment[totalnum.linesegment].end[0].gamma=q0.gamma;
	linesegment[totalnum.linesegment].end[1].v=q1.v;
	linesegment[totalnum.linesegment].end[1].w=q1.w;
	linesegment[totalnum.linesegment].end[1].gamma=q1.gamma;

	cresentsection[totalnum.cresentsection].numellipsearcs=1;
	cresentsection[totalnum.cresentsection].numscirclearcs=2;
	cresentsection[totalnum.cresentsection].linesegmentn=totalnum.linesegment;
	cresentsection[totalnum.cresentsection].ellipsearcn[0]=cresent[cren].ellipsearcn;
	cresentsection[totalnum.cresentsection].cw=(-1);
	scirclearc[totalnum.scirclearc].ellipsen=en;
	scirclearc[totalnum.scirclearc].radius=scirclearc[cresent[cren].scirclearcn].radius;
	scirclearc[totalnum.scirclearc].ebeta[0]=scirclearc[cresent[cren].scirclearcn].ebeta[0];
	scirclearc[totalnum.scirclearc].end[0].v=scirclearc[cresent[cren].scirclearcn].end[0].v;
	scirclearc[totalnum.scirclearc].end[0].w=scirclearc[cresent[cren].scirclearcn].end[0].w;
	scirclearc[totalnum.scirclearc].end[0].gamma=scirclearc[cresent[cren].scirclearcn].end[0].gamma;
	scirclearc[totalnum.scirclearc].ebeta[1]=circle[cirn].starintbeta[0];
	scirclearc[totalnum.scirclearc].end[1].v=q0.v;
	scirclearc[totalnum.scirclearc].end[1].w=q0.w;
	scirclearc[totalnum.scirclearc].end[1].gamma=q0.gamma;
	cresentsection[totalnum.cresentsection].scirclearcn[0]=totalnum.scirclearc;
	totalnum.scirclearc++;
	scirclearc[totalnum.scirclearc].ellipsen=en;
	scirclearc[totalnum.scirclearc].radius=scirclearc[cresent[cren].scirclearcn].radius;
	scirclearc[totalnum.scirclearc].ebeta[1]=scirclearc[cresent[cren].scirclearcn].ebeta[1];
	scirclearc[totalnum.scirclearc].end[1].v=scirclearc[cresent[cren].scirclearcn].end[1].v;
	scirclearc[totalnum.scirclearc].end[1].w=scirclearc[cresent[cren].scirclearcn].end[1].w;
	scirclearc[totalnum.scirclearc].end[1].gamma=scirclearc[cresent[cren].scirclearcn].end[1].gamma;
	scirclearc[totalnum.scirclearc].ebeta[0]=circle[cirn].starintbeta[1];
	scirclearc[totalnum.scirclearc].end[0].v=q1.v;
	scirclearc[totalnum.scirclearc].end[0].w=q1.w;
	scirclearc[totalnum.scirclearc].end[0].gamma=q1.gamma;
	cresentsection[totalnum.cresentsection].scirclearcn[1]=totalnum.scirclearc;
	totalnum.scirclearc++;
	cresentsection[totalnum.cresentsection].area=cresent[cren].area-areacirclesection2(starradius,circle[cirn].starintbeta[0],circle[cirn].starintbeta[1]);

	circlesection[totalnum.circlesection].linesegmentn=totalnum.linesegment;
	pcirclearc[totalnum.pcirclearc].circlen=cirn;
	pcirclearc[totalnum.pcirclearc].ellipsen=en;
	pcirclearc[totalnum.pcirclearc].end[0].v=q1.v;
	pcirclearc[totalnum.pcirclearc].end[0].w=q1.w;
	pcirclearc[totalnum.pcirclearc].end[0].gamma=q1.gamma;
	pcirclearc[totalnum.pcirclearc].end[1].v=q0.v;
	pcirclearc[totalnum.pcirclearc].end[1].w=q0.w;
	pcirclearc[totalnum.pcirclearc].end[1].gamma=q0.gamma;
	cwa=atan2(pcirclearc[totalnum.pcirclearc].end[0].w-wc,pcirclearc[totalnum.pcirclearc].end[0].v-vc);
	if(cwa<0)	cwa+=PIt2;
	ccwa=atan2(pcirclearc[totalnum.pcirclearc].end[1].w-wc,pcirclearc[totalnum.pcirclearc].end[1].v-vc);
	if(ccwa<0)	ccwa+=PIt2;
	if(ccwa<cwa)
		ccwa+=PIt2;
	pcirclearc[totalnum.pcirclearc].edang=ccwa-cwa;
	circlesection[totalnum.circlesection].pcirclearcn=totalnum.pcirclearc;
	totalnum.pcirclearc++;
	circlesection[totalnum.circlesection].area=areacirclesection(circle[cirn].radius,pcirclearc[circlesection[totalnum.circlesection].pcirclearcn].edang);

	totalnum.linesegment++;
	cresentpart[totalnum.cresentpart].circlesectionn=totalnum.circlesection;
	totalnum.circlesection++;
	cresentpart[totalnum.cresentpart].cresentsectionn=totalnum.cresentsection;
	totalnum.cresentsection++;
	cresentpart[totalnum.cresentpart].area=cresentsection[cresentpart[totalnum.cresentpart].cresentsectionn].area-circlesection[cresentpart[totalnum.cresentpart].circlesectionn].area;
	cresentpart[totalnum.cresentpart].ellipsen=en;
	totalnum.cresentpart++;
}
void createcresentpart11(vwpoint p,vwpoint q,char cw,int cren,int en,int cirn,double vc,double wc,double starradius)
{
	double ccwa,cwa;

	linesegment[totalnum.linesegment].ellipsen=en;
	if(cw)
	{
		linesegment[totalnum.linesegment].end[0].v=q.v;
		linesegment[totalnum.linesegment].end[0].w=q.w;
		linesegment[totalnum.linesegment].end[0].gamma=q.gamma;
		linesegment[totalnum.linesegment].end[1].v=p.v;
		linesegment[totalnum.linesegment].end[1].w=p.w;
		linesegment[totalnum.linesegment].end[1].gamma=p.gamma;
	}
	else
	{
		linesegment[totalnum.linesegment].end[0].v=p.v;
		linesegment[totalnum.linesegment].end[0].w=p.w;
		linesegment[totalnum.linesegment].end[0].gamma=p.gamma;
		linesegment[totalnum.linesegment].end[1].v=q.v;
		linesegment[totalnum.linesegment].end[1].w=q.w;
		linesegment[totalnum.linesegment].end[1].gamma=q.gamma;
	}

	cresentsection[totalnum.cresentsection].numellipsearcs=1;
	cresentsection[totalnum.cresentsection].numscirclearcs=1;
	cresentsection[totalnum.cresentsection].linesegmentn=totalnum.linesegment;
	cresentsection[totalnum.cresentsection].cw=cw;
	ellipsearc[totalnum.ellipsearc].ellipsen=en;
	if(cw)
	{
		ellipsearc[totalnum.ellipsearc].end[1].v=p.v;
		ellipsearc[totalnum.ellipsearc].end[1].w=p.w;
		ellipsearc[totalnum.ellipsearc].end[1].gamma=p.gamma;
		ellipsearc[totalnum.ellipsearc].end[0].v=ellipsearc[cresent[cren].ellipsearcn].end[0].v;
		ellipsearc[totalnum.ellipsearc].end[0].w=ellipsearc[cresent[cren].ellipsearcn].end[0].w;
		ellipsearc[totalnum.ellipsearc].end[0].gamma=ellipsearc[cresent[cren].ellipsearcn].end[0].gamma;
	}
	else
	{
		ellipsearc[totalnum.ellipsearc].end[0].v=p.v;
		ellipsearc[totalnum.ellipsearc].end[0].w=p.w;
		ellipsearc[totalnum.ellipsearc].end[0].gamma=p.gamma;
		ellipsearc[totalnum.ellipsearc].end[1].v=ellipsearc[cresent[cren].ellipsearcn].end[1].v;
		ellipsearc[totalnum.ellipsearc].end[1].w=ellipsearc[cresent[cren].ellipsearcn].end[1].w;
		ellipsearc[totalnum.ellipsearc].end[1].gamma=ellipsearc[cresent[cren].ellipsearcn].end[1].gamma;
	}
	cresentsection[totalnum.cresentsection].ellipsearcn[0]=totalnum.ellipsearc;
	totalnum.ellipsearc++;
	scirclearc[totalnum.scirclearc].ellipsen=en;
	scirclearc[totalnum.scirclearc].radius=scirclearc[cresent[cren].scirclearcn].radius;
	if(cw)
	{
		scirclearc[totalnum.scirclearc].ebeta[0]=scirclearc[cresent[cren].scirclearcn].ebeta[0];
		scirclearc[totalnum.scirclearc].end[0].v=scirclearc[cresent[cren].scirclearcn].end[0].v;
		scirclearc[totalnum.scirclearc].end[0].w=scirclearc[cresent[cren].scirclearcn].end[0].w;
		scirclearc[totalnum.scirclearc].end[0].gamma=scirclearc[cresent[cren].scirclearcn].end[0].gamma;
		scirclearc[totalnum.scirclearc].ebeta[1]=circle[cirn].starintbeta[0];
		scirclearc[totalnum.scirclearc].end[1].v=q.v;
		scirclearc[totalnum.scirclearc].end[1].w=q.w;
		scirclearc[totalnum.scirclearc].end[1].gamma=q.gamma;
	}
	else
	{
		scirclearc[totalnum.scirclearc].ebeta[0]=circle[cirn].starintbeta[1];
		scirclearc[totalnum.scirclearc].end[0].v=q.v;
		scirclearc[totalnum.scirclearc].end[0].w=q.w;
		scirclearc[totalnum.scirclearc].end[0].gamma=q.gamma;
		scirclearc[totalnum.scirclearc].ebeta[1]=scirclearc[cresent[cren].scirclearcn].ebeta[1];
		scirclearc[totalnum.scirclearc].end[1].v=scirclearc[cresent[cren].scirclearcn].end[1].v;
		scirclearc[totalnum.scirclearc].end[1].w=scirclearc[cresent[cren].scirclearcn].end[1].w;
		scirclearc[totalnum.scirclearc].end[1].gamma=scirclearc[cresent[cren].scirclearcn].end[1].gamma;
	}
	cresentsection[totalnum.cresentsection].scirclearcn[0]=totalnum.scirclearc;
	totalnum.scirclearc++;
	if(cw)
		cresentsection[totalnum.cresentsection].area=areatriangle(scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].end[1],ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[1],ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[0])
				+areacirclesection2(starradius,scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].ebeta[0],scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].ebeta[1])
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[0].gamma,ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[1].gamma);
	else
		cresentsection[totalnum.cresentsection].area=areatriangle(ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[0],scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].end[0],ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[1])
				+areacirclesection2(starradius,scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].ebeta[0],scirclearc[cresentsection[totalnum.cresentsection].scirclearcn[0]].ebeta[1])
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[0].gamma,ellipsearc[cresentsection[totalnum.cresentsection].ellipsearcn[0]].end[1].gamma);

	circlesection[totalnum.circlesection].linesegmentn=totalnum.linesegment;
	pcirclearc[totalnum.pcirclearc].circlen=cirn;
	pcirclearc[totalnum.pcirclearc].ellipsen=en;
	if(cw)
	{
		pcirclearc[totalnum.pcirclearc].end[0].v=p.v;
		pcirclearc[totalnum.pcirclearc].end[0].w=p.w;
		pcirclearc[totalnum.pcirclearc].end[0].gamma=p.gamma;
		pcirclearc[totalnum.pcirclearc].end[1].v=q.v;
		pcirclearc[totalnum.pcirclearc].end[1].w=q.w;
		pcirclearc[totalnum.pcirclearc].end[1].gamma=q.gamma;
	}
	else
	{
		pcirclearc[totalnum.pcirclearc].end[0].v=q.v;
		pcirclearc[totalnum.pcirclearc].end[0].w=q.w;
		pcirclearc[totalnum.pcirclearc].end[0].gamma=q.gamma;
		pcirclearc[totalnum.pcirclearc].end[1].v=p.v;
		pcirclearc[totalnum.pcirclearc].end[1].w=p.w;
		pcirclearc[totalnum.pcirclearc].end[1].gamma=p.gamma;
	}
	cwa=atan2(pcirclearc[totalnum.pcirclearc].end[0].w-wc,pcirclearc[totalnum.pcirclearc].end[0].v-vc);
	if(cwa<0)	cwa+=PIt2;
	ccwa=atan2(pcirclearc[totalnum.pcirclearc].end[1].w-wc,pcirclearc[totalnum.pcirclearc].end[1].v-vc);
	if(ccwa<0)	ccwa+=PIt2;
	if(ccwa<cwa)
		ccwa+=PIt2;
	pcirclearc[totalnum.pcirclearc].edang=ccwa-cwa;
	circlesection[totalnum.circlesection].pcirclearcn=totalnum.pcirclearc;
	totalnum.pcirclearc++;
	circlesection[totalnum.circlesection].area=areacirclesection(circle[cirn].radius,pcirclearc[circlesection[totalnum.circlesection].pcirclearcn].edang);

	totalnum.linesegment++;
	cresentpart[totalnum.cresentpart].circlesectionn=totalnum.circlesection;
	totalnum.circlesection++;
	cresentpart[totalnum.cresentpart].cresentsectionn=totalnum.cresentsection;
	totalnum.cresentsection++;
	cresentpart[totalnum.cresentpart].area=cresentsection[cresentpart[totalnum.cresentpart].cresentsectionn].area-circlesection[cresentpart[totalnum.cresentpart].circlesectionn].area;
	cresentpart[totalnum.cresentpart].ellipsen=en;
	if(cresentpart[totalnum.cresentpart].area<0)
	{	//this is a bit dodgy
		if(cw)	//endpoints of circle arc are wrong way round, try to fix them as best we can
		{
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].v=q.v;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].w=q.w;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].gamma=q.gamma;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].v=p.v;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].w=p.w;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].gamma=p.gamma;
		}
		else
		{
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].v=p.v;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].w=p.w;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[0].gamma=p.gamma;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].v=q.v;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].w=q.w;
			pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].end[1].gamma=q.gamma;
		}
		pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].edang=PIt2-pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].edang;
		circlesection[cresentpart[totalnum.cresentpart].circlesectionn].area=areacirclesection(circle[cirn].radius,pcirclearc[circlesection[cresentpart[totalnum.cresentpart].circlesectionn].pcirclearcn].edang);
		cresentpart[totalnum.cresentpart].area=cresentsection[cresentpart[totalnum.cresentpart].cresentsectionn].area-circlesection[cresentpart[totalnum.cresentpart].circlesectionn].area;
	}
	totalnum.cresentpart++;


}
char fixcresentnpnq(double starradius,int rn,int *np0,int *nq0,int wp[4],int wq[4])
{
	char pgood[4],qgood[4],minisp;
	int j,k,en,np=(*np0),nq=(*nq0),wmin;
	double dg0,dg1,pdang[4],qdang[4],tdang,mindang,ydif,zdif;
	vwpoint tp;

	en=cresent[rn].ellipsen;
	if(np>2)	//too many, try removing one
	{
		if(np==3)
		{	
			dg0=ellipse[en].circleint[0][0].gamma-PIo2;
			if(dg0<0) dg0=(-dg0);
			dg1=ellipse[en].circleint[0][1].gamma-PIo2;
			if(dg1<0) dg1=(-dg1);
			if(dg1>dg0)
			{
				dg0=dg1;
				k=1;
			}
			else
				k=0;
			dg1=ellipse[en].circleint[0][2].gamma-PIo2;
			if(dg1<0) dg1=(-dg1);
			if(dg1>dg0)
				k=2;
			for(j=k;j<2;j++)
				wp[j]=wp[j+1];
			np=2;
			(*np0)=2;
		}
		else
		{
			errorflag=4;
			return 0;
		}
	}

	if(np%2!=nq%2)
	{
		en=cresent[rn].ellipsen;
		for(j=0;j<ellipse[en].numcircleint[0];j++)	//find circle intersections with ellipse arc
		{
			if(isanglebetween(ellipsearc[cresent[rn].ellipsearcn].end[0].gamma,ellipsearc[cresent[rn].ellipsearcn].end[1].gamma,ellipse[en].circleint[0][j].gamma))
			{
				pgood[j]=1;
				pdang[j]=ellipse[en].circleint[0][j].gamma-ellipsearc[cresent[rn].ellipsearcn].end[0].gamma;
				if(pdang[j]<0) pdang[j]+=PIt2;
				tdang=ellipsearc[cresent[rn].ellipsearcn].end[1].gamma-ellipse[en].circleint[0][j].gamma;
				if(tdang<0) tdang+=PIt2;
				if(tdang<pdang[j])
					pdang[j]=tdang;
				if(wp[0]!=j&&wp[1]!=j)
					pgood[j]=0;
			}
			else
			{
				pgood[j]=0;
				pdang[j]=ellipsearc[cresent[rn].ellipsearcn].end[0].gamma-ellipse[en].circleint[0][j].gamma;
				if(pdang[j]<0) pdang[j]+=PIt2;
				tdang=ellipse[en].circleint[0][j].gamma-ellipsearc[cresent[rn].ellipsearcn].end[1].gamma;
				if(tdang<0) tdang+=PIt2;
				if(tdang<pdang[j])
					pdang[j]=tdang;
				if(wp[0]==j||wp[1]==j)
					return 0;
			}					
		}
		for(j=0;j<circle[0].numstarint;j++)	//find circle intersections with scircle arc
		{
			ydif=starradius*cos(circle[0].starintbeta[j])-ellipse[en].centery;
			zdif=starradius*sin(circle[0].starintbeta[j])-ellipse[en].centerz;
			tp.v=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	
			tp.w=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	
			tp.gamma=atan2(tp.w,tp.v);
			if(tp.gamma<0) tp.gamma+=PIt2;
			if(isanglebetween(scirclearc[cresent[rn].scirclearcn].ebeta[0],scirclearc[cresent[rn].scirclearcn].ebeta[1],circle[0].starintbeta[j]))
			{
				qgood[j]=1;
				qdang[j]=tp.gamma-scirclearc[cresent[rn].scirclearcn].end[0].gamma;
				if(qdang[j]<0) qdang[j]+=PIt2;
				tdang=scirclearc[cresent[rn].scirclearcn].end[1].gamma-tp.gamma;
				if(tdang<0) tdang+=PIt2;
				if(tdang<qdang[j])
					qdang[j]=tdang;
				if(wq[0]!=j&&wq[1]!=j)
					return 0;
			}
			else
			{
				qgood[j]=0;
				qdang[j]=scirclearc[cresent[rn].scirclearcn].end[0].gamma-tp.gamma;
				if(qdang[j]<0) qdang[j]+=PIt2;
				tdang=tp.gamma-scirclearc[cresent[rn].scirclearcn].end[1].gamma;
				if(tdang<0) tdang+=PIt2;
				if(tdang<qdang[j])
					qdang[j]=tdang;
				if(wq[0]==j||wq[1]==j)
					return 0;
			}
		}

		mindang=100.0;
		minisp=2;
		for(j=0;j<ellipse[en].numcircleint[0];j++)
			if(pdang[j]<mindang)
				if(pgood[j]==1||np<2)
				{
					mindang=pdang[j];
					wmin=j;
					minisp=1;
				}
		for(j=0;j<circle[0].numstarint;j++)
			if(qdang[j]<mindang)
				if(qgood[j]==1||nq<2)
				{
					mindang=qdang[j];
					wmin=j;
					minisp=0;
				}
		if(minisp==2)
			return 0;
		if(minisp)
		{
			if(pgood[wmin])
			{
				if(wp[0]==wmin)
					for(j=0;j<np-1;j++)
						wp[j]=wp[j+1];
				(*np0)--;
			}
			else
			{
				wp[np]=wmin;
				(*np0)++;
			}
		}
		else
		{
			if(qgood[wmin])
			{
				if(wq[0]=wmin)
					for(j=0;j<nq-1;j++)
						wq[j]=wq[j+1];
				(*nq0)--;
			}
			else
			{
				wq[nq]=wmin;
				(*nq0)++;
			}
		}
		if((*np0)%2!=(*nq0)%2)
			return 0;
		else
			return 1;
	}
	else
		return 1;

	return 0;
}
void checkoverlap(double starradius)
{	
	int i,j,k,np,nq,wp[4],wq[2],en;
	double ydif,zdif,dsq,dmax,vc,wc;
	vwpoint p[4],q[2];

	//for now, just one planet -> just one circle
	if(totalnum.circle>0)
	{
		for(i=0;i<totalnum.ellipse;i++)	//-----------check ellipses
		{
			ellipse[i].numcircleint[0]=0;
			ydif=circle[0].centery-ellipse[i].centery;
			zdif=circle[0].centerz-ellipse[i].centerz;
			dsq=ydif*ydif+zdif*zdif;
			dmax=circle[0].radius+ellipse[i].smajor;
			if(dsq<dmax*dmax)
			{
				vc=ydif*ellipse[i].vhaty+zdif*ellipse[i].vhatz;	//v-w coords of planet center
				wc=ydif*ellipse[i].whaty+zdif*ellipse[i].whatz;	//yes, really
				np=ellipsecircleintersection(ellipse[i].smajor,ellipse[i].sminor,vc,wc,circle[0].rsq,p);
				if(np==0)
				{
					if(ellipse[i].active)
					{
						if((vc-ellipse[i].smajor)*(vc-ellipse[i].smajor)+wc*wc<circle[0].rsq)	//point on ellipse inside circle
						{
							ellipse[i].active=0;	//ellipse completely covered
							whichspots[ellipse[i].spot]=1;
						}
						else if((vc+circle[0].radius)*(vc+circle[0].radius)/(ellipse[i].smajor*ellipse[i].smajor)+wc*wc/(ellipse[i].sminor*ellipse[i].sminor)<1.0)	//point on circle inside ellipse
						{	
							ellipsehole[totalnum.ellipsehole].ellipsen=i;
							ellipsehole[totalnum.ellipsehole].circlen=0;
							ellipsehole[totalnum.ellipsehole].area=ellipse[i].area-circle[0].area;
							totalnum.ellipsehole++;
							ellipse[i].active=0;
							whichspots[ellipse[i].spot]=1;
						}
					}
				}
				else if(np==2)
				{
					ellipse[i].numcircleint[0]=2;
					for(j=0;j<2;j++)
					{
						ellipse[i].circleint[0][j].v=p[j].v;	//ellipse under circle from 0 to 1
						ellipse[i].circleint[0][j].w=p[j].w;
						ellipse[i].circleint[0][j].gamma=p[j].gamma;	
					}
					if(ellipse[i].active)
					{
						createellipsepart(p[1],p[0],i,0,vc,wc);
						ellipse[i].active=0;
						whichspots[ellipse[i].spot]=1;
					}
				}
				else if(np==4)
				{
					ellipse[i].numcircleint[0]=4;
					for(j=0;j<4;j++)
					{
						ellipse[i].circleint[0][j].v=p[j].v;	//ellipse under circle: 0 to 1, 2 to 3
						ellipse[i].circleint[0][j].w=p[j].w;
						ellipse[i].circleint[0][j].gamma=p[j].gamma;
					}
					if(ellipse[i].active)
					{
						createellipsepart(p[1],p[2],i,0,vc,wc);
						createellipsepart(p[3],p[0],i,0,vc,wc);
						ellipse[i].active=0;
						whichspots[ellipse[i].spot]=1;
					}
				}
				else
					errorflag=4;
			}
		}


		for(i=0;i<totalnum.cresent;i++)	////-----------check cresents
		{
			k=0;
			en=cresent[i].ellipsen;
			np=0;	//find circle intersections with ellipse arc
			for(j=0;j<ellipse[en].numcircleint[0];j++)
				if(isanglebetween(ellipsearc[cresent[i].ellipsearcn].end[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma,ellipse[en].circleint[0][j].gamma))
				{
					wp[np]=j;
					np++;
				}
			nq=0;	//find circle intersections with scircle arc
			for(j=0;j<circle[0].numstarint;j++)
				if(isanglebetween(scirclearc[cresent[i].scirclearcn].ebeta[0],scirclearc[cresent[i].scirclearcn].ebeta[1],circle[0].starintbeta[j]))
				{
					if(nq>=2)
						errorflag=4;
					wq[nq]=j;
					nq++;
				}
			if(np%2!=nq%2||np>2)
			{
				if(!fixcresentnpnq(starradius,i,&np,&nq,wp,wq))
					errorflag=4;
			}

			ydif=circle[0].centery-ellipse[en].centery;
			zdif=circle[0].centerz-ellipse[en].centerz;
			vc=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	//v-w coords of planet center
			wc=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	//yes, really
			if(nq==0)
			{
				if(np==0)
				{
					if((ellipsearc[cresent[i].ellipsearcn].end[0].v-vc)*(ellipsearc[cresent[i].ellipsearcn].end[0].v-vc)+(ellipsearc[cresent[i].ellipsearcn].end[0].w-wc)*(ellipsearc[cresent[i].ellipsearcn].end[0].w-wc)<circle[0].rsq)
					{
						cresent[i].active=0;	//cresent covered by circle
						whichspots[ellipse[cresent[i].ellipsen].spot]=1;
					}
					else if(circle[0].mindcen>cresent[i].mindcen&&circle[0].maxdcen<scirclearc[cresent[i].scirclearcn].radius
							&&isanglebetween(ellipsearc[cresent[i].ellipsearcn].end[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma,atan2(wc,vc)))
					{	//circle is completly inside cresent
						cresenthole[totalnum.cresenthole].circlen=0;
						cresenthole[totalnum.cresenthole].cresentn=i;
						cresenthole[totalnum.cresenthole].area=cresent[i].area-circle[0].area;
						totalnum.cresenthole++;
						cresent[i].active=0;
						whichspots[ellipse[cresent[i].ellipsen].spot]=1;
					}
				}
				else if(np==2)	//np=2 nq=0
				{
					createcresentpart20(ellipse[en].circleint[0][wp[0]],ellipse[en].circleint[0][wp[1]],i,en,0,vc,wc);
					cresent[i].active=0;
					whichspots[ellipse[cresent[i].ellipsen].spot]=1;
				}
				else
					errorflag=4;
			}
			else if(nq==2)
			{
				for(j=0;j<2;j++)
				{
					ydif=starradius*cos(circle[0].starintbeta[wq[j]])-ellipse[en].centery;
					zdif=starradius*sin(circle[0].starintbeta[wq[j]])-ellipse[en].centerz;
					q[j].v=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	
					q[j].w=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	
					q[j].gamma=atan2(q[j].w,q[j].v);
					if(q[j].gamma<0) q[j].gamma+=PIt2;
				}
				if(np==0)	//np=0 nq=2
				{
					createcresentpart02(q[0],q[1],i,en,0,vc,wc,starradius);
					cresent[i].active=0;
					whichspots[ellipse[cresent[i].ellipsen].spot]=1;
				}
				else if(np==2)	//np=2 nq=2
				{
					createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],1,i,en,0,vc,wc,starradius);
					createcresentpart11(ellipse[en].circleint[0][wp[1]],q[1],0,i,en,0,vc,wc,starradius);
					cresent[i].active=0;
					whichspots[ellipse[cresent[i].ellipsen].spot]=1;
				}
				else
					errorflag=4;
			}
			else if(nq==1)
			{
				ydif=starradius*cos(circle[0].starintbeta[wq[0]])-ellipse[en].centery;
				zdif=starradius*sin(circle[0].starintbeta[wq[0]])-ellipse[en].centerz;
				q[0].v=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	
				q[0].w=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	
				q[0].gamma=atan2(q[0].w,q[0].v);
				if(q[0].gamma<0) q[0].gamma+=PIt2;
				if(np==1)
				{
					if(wq[0]%2!=wp[0]%2)
						errorflag=4;
					if(wq[0]==0)
						createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],1,i,en,0,vc,wc,starradius);
					else
						createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],0,i,en,0,vc,wc,starradius);
					cresent[i].active=0;
					whichspots[ellipse[cresent[i].ellipsen].spot]=1;
				}
				else
					errorflag=4;
			}
			else
				errorflag=4;
		}
	}
}
void checkstarcircleint(double starradius,int cn)
{
	double starradiussq,a,c,csq,d,anga,angb;

	if(starradius<circle[cn].mindcen||starradius>circle[cn].maxdcen)
		circle[cn].numstarint=0;
	else
	{
		d=sqrt(circle[cn].centery*circle[cn].centery+circle[cn].centerz*circle[cn].centerz);
		starradiussq=starradius*starradius;
		a=(starradiussq-circle[cn].rsq+d*d)/(2*d);
		csq=starradiussq-a*a;
		if(csq<0)
			circle[cn].numstarint=0;
		else
		{
			c=sqrt(csq);
			angb=atan2(circle[cn].centerz,circle[cn].centery);
			anga=acos(a/starradius);
			circle[cn].starintbeta[0]=angb-anga;
			if(circle[cn].starintbeta[0]<0)	circle[cn].starintbeta[0]+=PIt2;
			circle[cn].starintbeta[1]=angb+anga;
			if(circle[cn].starintbeta[1]<0)	circle[cn].starintbeta[1]+=PIt2;
			a=asin(c/circle[cn].radius);
			if(d*d<starradiussq-circle[cn].rsq)
				a=PI-a;	//mostly inside real star
			circle[cn].starintdang=2.0*a;
			circle[cn].scy[0]=starradius*cos(circle[cn].starintbeta[0]);
			circle[cn].scz[0]=starradius*sin(circle[cn].starintbeta[0]);
			circle[cn].scy[1]=starradius*cos(circle[cn].starintbeta[1]);
			circle[cn].scz[1]=starradius*sin(circle[cn].starintbeta[1]);
			circle[cn].numstarint=2;
		}
	}
}
double realstararea(double starradius)
{
	int i;
	double totarea,darea;

#	if PRINTEACHSPOT
		for(i=0;i<numspots;i++)
			ldspotreport[i]=0.0;
#	endif
	totarea=PI*starradius*starradius;
	for(i=0;i<totalnum.circle;i++)
	{
		checkstarcircleint(starradius,i);
		if(circle[i].numstarint==0)
		{
			if(starradius>circle[i].mindcen)
				if(starradius<circle[i].maxdcen)
					totarea-=PI*starradius*starradius;
				else
					totarea-=circle[i].area;	
		}
		else
		{
			totarea-=areacirclesection(circle[i].radius,circle[i].starintdang);
			totarea-=areacirclesection2(starradius,circle[i].starintbeta[0],circle[i].starintbeta[1]);
		}
	}
	for(i=0;i<totalnum.ellipse;i++)
		if(ellipse[i].active)
		{
			darea=ellipse[i].area*ellipse[i].fracdarkness;
#			if PRINTEACHSPOT
				ldspotreport[ellipse[i].spot]+=darea;
#			endif
			totarea-=darea;
		}
	for(i=0;i<totalnum.cresent;i++)
		if(cresent[i].active)
		{
			darea=cresent[i].area*ellipse[cresent[i].ellipsen].fracdarkness;
#			if PRINTEACHSPOT 
				ldspotreport[ellipse[cresent[i].ellipsen].spot]+=darea;
#			endif
			totarea-=darea;
		}
	for(i=0;i<totalnum.ellipsepart;i++)
	{
		darea=ellipsepart[i].area*ellipse[ellipsepart[i].ellipsen].fracdarkness;
#		if PRINTEACHSPOT 
			ldspotreport[ellipse[ellipsepart[i].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}
	for(i=0;i<totalnum.cresentpart;i++)
	{
			darea=cresentpart[i].area*ellipse[cresentpart[i].ellipsen].fracdarkness;
#			if PRINTEACHSPOT
				ldspotreport[ellipse[cresentpart[i].ellipsen].spot]+=darea;
#			endif
			totarea-=darea;
	}
	for(i=0;i<totalnum.ellipsehole;i++)
	{
		darea=ellipsehole[i].area*ellipse[ellipsehole[i].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[ellipsehole[i].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}
	for(i=0;i<totalnum.cresenthole;i++)
	{
		darea=cresenthole[i].area*ellipse[cresent[cresenthole[i].cresentn].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[cresent[cresenthole[i].cresentn].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}
	if(totarea<0||totarea>PI*starradius*starradius)
		errorflag=5;
	return totarea;
}
double ldcirclearea(double starradius,int cn)
{
	if(circle[cn].numstarint==0)
	{
		if(starradius<circle[cn].mindcen)
			return 0.0;	
		else if(starradius<circle[cn].maxdcen)
			return PI*starradius*starradius;
		else
			return circle[cn].area;	
	}
	return areacirclesection(circle[cn].radius,circle[cn].starintdang)+areacirclesection2(starradius,circle[cn].starintbeta[0],circle[cn].starintbeta[1]);
}
void checkstarellipseint(double starradius,int en)
{
	int n;
	double v[2],w[2],gamma[2];

	if(starradius>ellipse[en].mindcen&&starradius<ellipse[en].maxdcen)
	{
		n=ellipsestarintersection(starradius,ellipse[en].dstarcenter,ellipse[en].smajor,ellipse[en].sminor,v,w);
		if(n==0)
			ellipse[en].numstarint=0;
		else if(n==1)	//1 pair
		{
			ellipse[en].numstarint=2;
			gamma[0]=atan2(w[0],v[0]);
			if(gamma[0]<0)	gamma[0]+=PIt2;
			if(ellipse[en].dstarcenter+ellipse[en].sminor>starradius)
			{
				ellipse[en].se[0].v=v[0];
				ellipse[en].se[0].w=w[0];
				ellipse[en].se[0].gamma=gamma[0];
				ellipse[en].se[1].v=(-v[0]);
				ellipse[en].se[1].w=w[0];
				ellipse[en].se[1].gamma=PI-gamma[0];
				if(ellipse[en].se[1].gamma<0)	ellipse[en].se[1].gamma+=PIt2;
			}
			else
			{
				ellipse[en].se[1].v=v[0];
				ellipse[en].se[1].w=w[0];
				ellipse[en].se[1].gamma=gamma[0];
				ellipse[en].se[0].v=(-v[0]);
				ellipse[en].se[0].w=w[0];
				ellipse[en].se[0].gamma=PI-gamma[0];
				if(ellipse[en].se[0].gamma<0)	ellipse[en].se[0].gamma+=PIt2;
			}
		}
		else if(n==2)	//2 pairs
		{
			ellipse[en].numstarint=4;
			gamma[0]=atan2(w[0],v[0]);
			if(gamma[0]<0)	gamma[0]+=PIt2;
			gamma[1]=atan2(w[1],v[1]);
			if(gamma[1]<0)	gamma[1]+=PIt2;
			ellipse[en].se[0].v=v[0];
			ellipse[en].se[0].w=w[0];
			ellipse[en].se[0].gamma=gamma[0];
			ellipse[en].se[1].v=v[1];
			ellipse[en].se[1].w=w[1];
			ellipse[en].se[1].gamma=gamma[1];
			ellipse[en].se[3].v=(-v[0]);
			ellipse[en].se[3].w=w[0];
			ellipse[en].se[3].gamma=PI-gamma[0];
			if(ellipse[en].se[3].gamma<0)	ellipse[en].se[3].gamma+=PIt2;
			ellipse[en].se[2].v=(-v[1]);
			ellipse[en].se[2].w=w[1];
			ellipse[en].se[2].gamma=PI-gamma[1];
			if(ellipse[en].se[2].gamma<0)	ellipse[en].se[2].gamma+=PIt2;
		}
		else
			errorflag=6;
	}
	else
		ellipse[en].numstarint=0;
}
double ldellipsearea(double starradius,int en)
{
	if(ellipse[en].numstarint==0)
	{
		if(starradius<ellipse[en].dstarcenter-ellipse[en].sminor)
			return 0.0;	//ellipse and star are disjoint
		if(starradius-ellipse[en].dstarcenter>ellipse[en].sminor)
			return ellipse[en].area;	//ellipse is in star
		return PI*starradius*starradius;	//star is in ellipse
	}
	if(ellipse[en].numstarint==2)
		return areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[1].gamma,ellipse[en].se[0].gamma)
				+areacirclesection2(starradius,sebeta(0,en),sebeta(1,en));
	if(ellipse[en].numstarint==4)
		return ellipse[en].area-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[0].gamma,ellipse[en].se[1].gamma)
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[2].gamma,ellipse[en].se[3].gamma)
				+areacirclesection2(starradius,sebeta(0,en),sebeta(1,en))
				+areacirclesection2(starradius,sebeta(2,en),sebeta(3,en));
	errorflag=7;
	return 0.0;
}
double ldcresentarea(double starradius,int rn)
{
	int i,nsr,sr[2],en;
	
	en=cresent[rn].ellipsen;
	nsr=0;
	if(ellipse[en].numstarint==2)
	{
		for(i=1;i>=0;i--)	//count backwards to get order right
			if(isanglebetween(ellipsearc[cresent[rn].ellipsearcn].end[0].gamma,ellipsearc[cresent[rn].ellipsearcn].end[1].gamma,ellipse[en].se[i].gamma))
			{
				if(nsr>=2)
				{
					errorflag=8;
					return 0.0;
				}
				sr[nsr]=i;	
				nsr++;
			}
	}
	else if(ellipse[en].numstarint==4)
	{
		for(i=0;i<4;i++)
			if(isanglebetween(ellipsearc[cresent[rn].ellipsearcn].end[0].gamma,ellipsearc[cresent[rn].ellipsearcn].end[1].gamma,ellipse[en].se[i].gamma))
			{
				if(nsr>=2)
				{
					errorflag=8;
					return 0.0;
				}
				sr[nsr]=i;	
				nsr++;
			}
	}

	if(nsr==0)
		return 0.0;
	if(nsr==2)
		return areacirclesection2(starradius,sebeta(sr[0],en),sebeta(sr[1],en))
			-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[sr[0]].gamma,ellipse[en].se[sr[1]].gamma);
	if(nsr==1)
	{
		fprintf(outerr,"one cresent intersection\n");
		errorflag=108;
		return 0.0;
	}
	errorflag=8;
	return 0.0;
}
int checkstarlinesegmentint(double starradius,int ln,vwpoint sl[2])
{	//finds star-segment intersections, sl[0] is closer to end[0]
	int nsl;
	double sv,sw,p0wpd,a,b,c,disc,tp,tn;

	sv=linesegment[ln].end[1].v-linesegment[ln].end[0].v;
	sw=linesegment[ln].end[1].w-linesegment[ln].end[0].w;
	p0wpd=linesegment[ln].end[0].w+ellipse[linesegment[ln].ellipsen].dstarcenter;
	a=sv*sv+sw*sw;
	b=2.0*(linesegment[ln].end[0].v*sv+(p0wpd)*sw);
	c=linesegment[ln].end[0].v*linesegment[ln].end[0].v+p0wpd*p0wpd-starradius*starradius;
	disc=b*b-4.0*a*c;
	if(disc<0)
		return 0;
	disc=sqrt(disc);
	a*=2.0;
	tp=(-b+disc)/a;
	tn=(-b-disc)/a;
	if(tn>=0.0&&tn<=1.0)
	{
		sl[0].v=linesegment[ln].end[0].v+tn*sv;
		sl[0].w=linesegment[ln].end[0].w+tn*sw;
		nsl=1;
	}
	else
		nsl=0;
	if(tp>=0.0&&tp<=1.0&&tp!=tn)
	{
		sl[nsl].v=linesegment[ln].end[0].v+tp*sv;
		sl[nsl].w=linesegment[ln].end[0].w+tp*sw;
		nsl++;
	}
	return nsl;
}
double ldellipsesectionarea(double starradius,int esn,vwpoint sl[2],int nsl)
{
	int i,en,wse[4],nse;
	double ydif,zdif,wdif;
	double starradiussq;

	starradiussq=starradius*starradius;
	en=ellipsearc[ellipsesection[esn].ellipsearcn].ellipsen;

	if(ellipse[en].numstarint==0)
	{
		if(starradius<ellipse[en].dstarcenter-ellipse[en].sminor)
			return 0.0;	//ellipse and star are disjoint
		if(starradius-ellipse[en].dstarcenter>ellipse[en].sminor)
			return ellipsesection[esn].area;	//ellipse section is in star
		nse=0;
	}
	else if(ellipse[en].numstarint==2)
	{
		nse=0;
		for(i=0;i<2;i++)
			if(isanglebetween(ellipsearc[ellipsesection[esn].ellipsearcn].end[0].gamma,ellipsearc[ellipsesection[esn].ellipsearcn].end[1].gamma,ellipse[en].se[i].gamma))
			{
				wse[nse]=i;
				nse++;
			}
	}
	else
	{
		nse=0;
		for(i=0;i<4;i++)
			if(isanglebetween(ellipsearc[ellipsesection[esn].ellipsearcn].end[0].gamma,ellipsearc[ellipsesection[esn].ellipsearcn].end[1].gamma,ellipse[en].se[i].gamma))
			{
				wse[nse]=i;
				nse++;
			}
		if(nse==2&&wse[0]==1&&wse[1]==2)
		{
			wse[0]=2;
			wse[1]=1;
		}
		else if(nse==3&&wse[1]==2&&wse[2]==3)
		{
			wse[2]=wse[0];
			wse[1]=3;
			wse[0]=2;
		}
	}

	if(nsl==0)
	{
		if(nse==0)
		{
			vwtoyz(ellipsearc[ellipsesection[esn].ellipsearcn].end[0].v,ellipsearc[ellipsesection[esn].ellipsearcn].end[0].w,en,&ydif,&zdif);
			if(ydif*ydif+zdif*zdif<starradiussq)	//point on ellipse section is inside star
				return ellipsesection[esn].area;	//all ellipse section is in star
			wdif=starradius-ellipse[en].dstarcenter;
			if(wdif>(-ellipse[en].sminor)&&wdif<ellipse[en].sminor)	//point on star edge is inside ellipse
				if(wonleft(ellipsesection[esn].linesegmentn,wdif))	//point on star is left of the line
					return PI*starradiussq;	//all star is inside ellipse section
			return 0.0;	//they are disjoint
		}
		else if(nse==2)
		{
			if(isanglebetween(ellipsearc[ellipsesection[esn].ellipsearcn].end[0].gamma,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma))
				return areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[1]].gamma,ellipse[en].se[wse[0]].gamma)
					+areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en));
			else
				return ellipsesection[esn].area
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma)
					+areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en));
		}
		else if(nse==4)	//yes, original ellipse area
			return ellipse[en].area-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[0].gamma,ellipse[en].se[1].gamma)
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[2].gamma,ellipse[en].se[3].gamma)
				+areacirclesection2(starradius,sebeta(0,en),sebeta(1,en))
				+areacirclesection2(starradius,sebeta(2,en),sebeta(3,en));
		else
		{
			errorflag=110;
			return 0.0;
		}
	}
	else if(nsl==1)
	{
		vwtoyz(linesegment[ellipsesection[esn].linesegmentn].end[0].v,linesegment[ellipsesection[esn].linesegmentn].end[0].w,en,&ydif,&zdif);
		if(nse==1)
		{
			if(ydif*ydif+zdif*zdif<starradiussq)	//end 0 is in star
				return areatriangle(linesegment[ellipsesection[esn].linesegmentn].end[0],sl[0],ellipse[en].se[wse[0]])
					+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,linesegment[ellipsesection[esn].linesegmentn].end[0].gamma)
					+areacirclesection2(starradius,vwpointtobeta(sl[0],en),sebeta(wse[0],en));
			else
			{
				vwtoyz(linesegment[ellipsesection[esn].linesegmentn].end[1].v,linesegment[ellipsesection[esn].linesegmentn].end[1].w,en,&ydif,&zdif);
				if(ydif*ydif+zdif*zdif<starradiussq)	//end 1 is in star
					return areatriangle(sl[0],linesegment[ellipsesection[esn].linesegmentn].end[1],ellipse[en].se[wse[0]])
						+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[ellipsesection[esn].linesegmentn].end[1].gamma,ellipse[en].se[wse[0]].gamma)
						+areacirclesection2(starradius,sebeta(wse[0],en),vwpointtobeta(sl[0],en));
				else
				{
					errorflag=110;
					return 0.0;
				}
			}
		}
		else if(nse==3)
		{
			if(ydif*ydif+zdif*zdif<starradiussq)	//end 0 is in star
			{
				return ellipsesection[esn].area-areatriangle(sl[0],linesegment[ellipsesection[esn].linesegmentn].end[1],ellipse[en].se[wse[2]])
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[ellipsesection[esn].linesegmentn].end[1].gamma,ellipse[en].se[wse[2]].gamma)
					+areacirclesection2(starradius,vwpointtobeta(sl[0],en),sebeta(wse[2],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma)
					+areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en));
			}
			else	//end 1 is in star
			{
				return ellipsesection[esn].area-areatriangle(linesegment[ellipsesection[esn].linesegmentn].end[0],sl[0],ellipse[en].se[wse[2]])
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[2]].gamma,linesegment[ellipsesection[esn].linesegmentn].end[0].gamma)
					+areacirclesection2(starradius,sebeta(wse[2],en),vwpointtobeta(sl[0],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma)
					+areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en));
			}
		}
		else
		{
			errorflag=110;
			return 0.0;
		}
	}
	else if(nsl==2)
	{
		if(nse==0)
			return areacirclesection2(starradius,vwpointtobeta(sl[1],en),vwpointtobeta(sl[0],en));
		else if(nse==2)
			return ellipsesection[esn].area-areatriangle(linesegment[ellipsesection[esn].linesegmentn].end[0],sl[0],ellipse[en].se[wse[0]])
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,linesegment[ellipsesection[esn].linesegmentn].end[0].gamma)
				+areacirclesection2(starradius,sebeta(wse[0],en),vwpointtobeta(sl[0],en))
				-areatriangle(sl[1],linesegment[ellipsesection[esn].linesegmentn].end[1],ellipse[en].se[wse[1]])
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[ellipsesection[esn].linesegmentn].end[1].gamma,ellipse[en].se[wse[1]].gamma)
				+areacirclesection2(starradius,vwpointtobeta(sl[1],en),sebeta(wse[1],en));
		else if(nse==4)
			return areaquad(ellipse[en].se)
			+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[1]].gamma,ellipse[en].se[wse[2]].gamma)
			+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[3]].gamma,ellipse[en].se[wse[0]].gamma)
			+areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en))
			+areacirclesection2(starradius,sebeta(wse[2],en),sebeta(wse[3],en))
			-areacirclesection2(starradius,vwpointtobeta(sl[0],en),vwpointtobeta(sl[1],en));
		else
		{
			errorflag=110;
			return 0.0;
		}
	}
	errorflag=10;
	return 0.0;
}
double ldcresentsectionarea(double starradius,int rsn,vwpoint sl[2],int nsl)
{
	int i,j,en,wa[4],wse[4],nsr;
	double ydif,zdif;
	double starradiussq;

	starradiussq=starradius*starradius;
	en=ellipsearc[cresentsection[rsn].ellipsearcn[0]].ellipsen;
	if(ellipse[en].numstarint==0)
	{
		if(starradius<ellipse[en].dstarcenter-ellipse[en].sminor)
			return 0.0;	//ellipse and star are disjoint
		nsr=0;
	}
	else if(ellipse[en].numstarint==2)
	{
		nsr=0;
		for(i=1;i>=0;i--)	//count backwards so star is in cresent from 0 to 1
			for(j=0;j<cresentsection[rsn].numellipsearcs;j++)
				if(isanglebetween(ellipsearc[cresentsection[rsn].ellipsearcn[j]].end[0].gamma,ellipsearc[cresentsection[rsn].ellipsearcn[j]].end[1].gamma,ellipse[en].se[i].gamma))
				{
					wse[nsr]=i;
					wa[nsr]=j;
					nsr++;
				}

	}
	else
	{
		nsr=0;
		for(i=0;i<4;i++)	//if a pair, should be se[1],se[2]
			for(j=0;j<cresentsection[rsn].numellipsearcs;j++)
				if(isanglebetween(ellipsearc[cresentsection[rsn].ellipsearcn[j]].end[0].gamma,ellipsearc[cresentsection[rsn].ellipsearcn[j]].end[1].gamma,ellipse[en].se[i].gamma))
				{
					wse[nsr]=i;
					wa[nsr]=j;
					nsr++;
				}
	}

	if(nsl==0)
	{
		if(nsr==0)
			return 0.0;
		else if(nsr==2)
		{
			if(wa[0]==wa[1])
				return areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma);
			else
				return areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma)
					+areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[cresentsection[rsn].linesegmentn].end[1].gamma,linesegment[cresentsection[rsn].linesegmentn].end[0].gamma);
		}
		else
		{
			if(cresentsection[rsn].area>ACCEPTABLEERROR)
				errorflag=112;
			return 0.0;
		}
	}
	else if(nsl==1)
	{
		if(nsr!=1)
		{
			if(cresentsection[rsn].area>ACCEPTABLEERROR)
				errorflag=112;
			return 0.0;
		}
		vwtoyz(linesegment[cresentsection[rsn].linesegmentn].end[0].v,linesegment[cresentsection[rsn].linesegmentn].end[0].w,en,&ydif,&zdif);
		if(ydif*ydif+zdif*zdif<starradiussq)	//end 0 is in star
		{
			if(cresentsection[rsn].numellipsearcs==1)
				return areatriangle(linesegment[cresentsection[rsn].linesegmentn].end[0],sl[0],ellipse[en].se[wse[0]])
					+areacirclesection2(starradius,vwpointtobeta(sl[0],en),sebeta(wse[0],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[cresentsection[rsn].linesegmentn].end[0].gamma,ellipse[en].se[wse[0]].gamma);
			else
				return areacirclesection2(starradius,vwpointtobeta(sl[0],en),sebeta(wse[0],en))
					-areatriangle(sl[0],linesegment[cresentsection[rsn].linesegmentn].end[0],ellipse[en].se[wse[0]])
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,linesegment[cresentsection[rsn].linesegmentn].end[0].gamma,ellipse[en].se[wse[0]].gamma);
		}
		else	//end 1 is in star
		{
			if(cresentsection[rsn].numellipsearcs==1)
				return areatriangle(sl[0],linesegment[cresentsection[rsn].linesegmentn].end[1],ellipse[en].se[wse[0]])
					+areacirclesection2(starradius,sebeta(wse[0],en),vwpointtobeta(sl[0],en))
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,linesegment[cresentsection[rsn].linesegmentn].end[1].gamma);
			else
				return areacirclesection2(starradius,sebeta(wse[0],en),vwpointtobeta(sl[0],en))
					-areatriangle(linesegment[cresentsection[rsn].linesegmentn].end[1],sl[0],ellipse[en].se[wse[0]])
					-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,linesegment[cresentsection[rsn].linesegmentn].end[1].gamma);
		}
	}
	else if(nsl==2)
	{
		if(nsr==2&&cresentsection[rsn].numellipsearcs==1)
			return areacirclesection2(starradius,sebeta(wse[0],en),sebeta(wse[1],en))
				-areacirclesection2(starradius,vwpointtobeta(sl[0],en),vwpointtobeta(sl[1],en))
				-areaellipsesection(ellipse[en].smajor,ellipse[en].sminor,ellipse[en].se[wse[0]].gamma,ellipse[en].se[wse[1]].gamma);
		else if(nsr==0&&cresentsection[rsn].numellipsearcs==2)
			return areacirclesection2(starradius,vwpointtobeta(sl[1],en),vwpointtobeta(sl[0],en));
		else
		{
			if(cresentsection[rsn].area>ACCEPTABLEERROR)
				errorflag=112;
			return 0.0;
		}
	}
	if(cresentsection[rsn].area>ACCEPTABLEERROR)
		errorflag=112;
	return 0.0;
}
double ldcirclesectionarea(double starradius,int csn,vwpoint sl[2],int nsl)
{
	int i,cn,nsc;
	double ta,ty,tz,l1y,l1z,starradiussq,scbeta[2];
	vwpoint sc[2],slloc[2],l[2];	//these are in the circle reference frame (vhat=yhat,what=zhat,centered on circle)

	starradiussq=starradius*starradius;
	cn=pcirclearc[circlesection[csn].pcirclearcn].circlen;
	if(circle[cn].numstarint==0)
	{
		if(starradius<circle[cn].mindcen)
			return 0.0;	
		else if(starradius>circle[cn].maxdcen)
			return circlesection[csn].area;
		nsc=0;
	}
	else if(circle[cn].numstarint==2)
	{
		for(i=0;i<2;i++)
		{
			vwtoyz(pcirclearc[circlesection[csn].pcirclearcn].end[1-i].v,pcirclearc[circlesection[csn].pcirclearcn].end[1-i].w,linesegment[circlesection[csn].linesegmentn].ellipsen,&ty,&tz);
			l1y=ty;	//for later
			l1z=tz;	//
			l[i].v=ty-circle[cn].centery;
			l[i].w=tz-circle[cn].centerz;
			l[i].gamma=atan2(l[i].w,l[i].v);
			if(l[i].gamma<0) l[i].gamma+=PIt2;
			sc[i].v=circle[cn].scy[i]-circle[cn].centery;
			sc[i].w=circle[cn].scz[i]-circle[cn].centerz;
			sc[i].gamma=atan2(sc[i].w,sc[i].v);
			if(sc[i].gamma<0) sc[i].gamma+=PIt2;
			scbeta[i]=circle[cn].starintbeta[i];
		}
		if(isanglebetween(l[1].gamma,l[0].gamma,sc[0].gamma))	//note l's are numbered backwards from arc ends
			nsc=1;
		else
			nsc=0;
		if(nsc==1&&isanglebetween(l[1].gamma,l[0].gamma,sc[1].gamma))
			nsc=2;
		else if(nsc==0&&isanglebetween(l[1].gamma,l[0].gamma,sc[1].gamma))
		{
			nsc=1;
			sc[0].v=sc[1].v; sc[0].w=sc[1].w; sc[0].gamma=sc[1].gamma; scbeta[0]=scbeta[1];
		}
	}
	else
	{
		errorflag=115;
		return 0.0;
	}

	for(i=0;i<nsl;i++)
	{
		vwtoyz(sl[i].v,sl[i].w,linesegment[circlesection[csn].linesegmentn].ellipsen,&ty,&tz);
		slloc[i].v=ty-circle[cn].centery;
		slloc[i].w=tz-circle[cn].centerz;
		slloc[i].gamma=atan2(slloc[i].w,slloc[i].v);
		if(slloc[i].gamma<0) slloc[i].gamma+=PIt2;
	}

	if(nsl==0)
	{
		if(nsc==0)
		{
			vwtoyz(linesegment[circlesection[csn].linesegmentn].end[0].v,linesegment[circlesection[csn].linesegmentn].end[0].w,linesegment[circlesection[csn].linesegmentn].ellipsen,&ty,&tz);	
			if(ty*ty+tz*tz<starradiussq)	//point on section is in star
				return circlesection[csn].area;	//whole section is in star
			if(starradius>circle[cn].mindcen&&starradius<circle[cn].maxdcen)	//point on star is inside circle
			{
				if(sqrt(circle[cn].centery*circle[cn].centery+circle[cn].centerz*circle[cn].centerz)+starradius>circle[cn].radius)
					return 0.0; //disjoint
				tz=starradius-ellipse[linesegment[circlesection[csn].linesegmentn].ellipsen].dstarcenter;
				if(wonleft(circlesection[csn].linesegmentn,tz))	//and on left of line segment
					return PI*starradiussq;
			}
			return 0.0;	//they are disjoint
		}
		else if(nsc==2)
			return areacirclesection2(circle[cn].radius,sc[1].gamma,sc[0].gamma)+areacirclesection2(starradius,scbeta[0],scbeta[1]);
		else
		{
			errorflag=115;
			return 0.0;
		}
	}
	else if(nsl==1)
	{
		if(nsc!=1)
		{
			errorflag=115;
			return 0.0;
		}
		if(l1y*l1y+l1z*l1z<starradiussq)	//end 1 of segment is inside star
		{
			ta=areatriangle(slloc[0],l[1],sc[0]);
			if(ta>=0)
				return ta+areacirclesection2(starradius,scbeta[0],vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen))
					+areacirclesection2(circle[cn].radius,l[1].gamma,sc[0].gamma);
			else	//circle points are in wrong order
				return (-ta)+areacirclesection2(starradius,vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen),scbeta[0])
					+areacirclesection2(circle[cn].radius,l[1].gamma,sc[0].gamma);
		}
		else	//end 0 of segment is in star
		{
			ta=areatriangle(l[0],slloc[0],sc[0]);
			if(ta>=0)
				return ta+areacirclesection2(starradius,vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen),scbeta[0])
					+areacirclesection2(circle[cn].radius,sc[0].gamma,l[0].gamma);
			else	//circle points in wrong order
				return (-ta)+areacirclesection2(starradius,scbeta[0],vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen))
					+areacirclesection2(circle[cn].radius,sc[0].gamma,l[0].gamma);
		}
	}
	else if(nsl==2)
	{
		if(nsc==0)
			return areacirclesection2(starradius,vwpointtobeta(sl[1],linesegment[circlesection[csn].linesegmentn].ellipsen),vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen));
		else if(nsc==2)
			return areacirclesection2(starradius,vwpointtobeta(sl[1],linesegment[circlesection[csn].linesegmentn].ellipsen),vwpointtobeta(sl[0],linesegment[circlesection[csn].linesegmentn].ellipsen))
				-areacirclesection2(starradius,scbeta[1],scbeta[0])
				+areacirclesection2(circle[cn].radius,sc[1].gamma,sc[0].gamma);
		else
		{
			errorflag=115;
			return 0.0;
		}
	}
	errorflag=115;
	return 0.0;
}
double ldellipsepartarea(double starradius,int epn)
{
	int nsl;
	double esarea;
	vwpoint sl[2];

	nsl=checkstarlinesegmentint(starradius,ellipsesection[ellipsepart[epn].ellipsesectionn].linesegmentn,sl);
	
 	esarea=ldellipsesectionarea(starradius,ellipsepart[epn].ellipsesectionn,sl,nsl);
	if(esarea==0.0)
		return 0.0;
	else
		return esarea-ldcirclesectionarea(starradius,ellipsepart[epn].circlesectionn,sl,nsl);
}
double ldcresentpartarea(double starradius,int rpn)
{
	int nsl;
	double csarea;
	vwpoint sl[2];

	nsl=checkstarlinesegmentint(starradius,cresentsection[cresentpart[rpn].cresentsectionn].linesegmentn,sl);

	csarea=ldcresentsectionarea(starradius,cresentpart[rpn].cresentsectionn,sl,nsl);
	if(csarea==0.0)
		return 0.0;
	else
		return csarea-ldcirclesectionarea(starradius,cresentpart[rpn].circlesectionn,sl,nsl);
}
double ldstararea(double starradius)
{
	int i;
	double totarea,darea;

#	if PRINTEACHSPOT
		for(i=0;i<numspots;i++)
			ldspotreport[i]=0.0;
#	endif
	totarea=PI*starradius*starradius;

	for(i=0;i<totalnum.circle;i++)
	{
		checkstarcircleint(starradius,i);
		totarea-=ldcirclearea(starradius,i);
	}

	for(i=0;i<totalnum.ellipse;i++)
	{
		checkstarellipseint(starradius,i);
		if(ellipse[i].active)
		{
			darea=ldellipsearea(starradius,i)*ellipse[i].fracdarkness;
#			if PRINTEACHSPOT
				ldspotreport[ellipse[i].spot]+=darea;
#			endif
			totarea-=darea;
		}
	}

	for(i=0;i<totalnum.cresent;i++)
	{
		if(cresent[i].active)
		{
			darea=ldcresentarea(starradius,i)*ellipse[cresent[i].ellipsen].fracdarkness;
#			if PRINTEACHSPOT
				ldspotreport[ellipse[cresent[i].ellipsen].spot]+=darea;
#			endif
			totarea-=darea;
		}
	}

	for(i=0;i<totalnum.ellipsehole;i++)
	{
		darea=(ldellipsearea(starradius,ellipsehole[i].ellipsen)-ldcirclearea(starradius,ellipsehole[i].circlen))*ellipse[ellipsehole[i].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[ellipsehole[i].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}

	for(i=0;i<totalnum.cresenthole;i++)
	{
		darea=(ldcresentarea(starradius,cresenthole[i].cresentn)-ldcirclearea(starradius,cresenthole[i].circlen))*ellipse[cresent[cresenthole[i].cresentn].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[cresent[cresenthole[i].cresentn].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}

	for(i=0;i<totalnum.ellipsepart;i++)
	{
		darea=ldellipsepartarea(starradius,i)*ellipse[ellipsepart[i].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[ellipsepart[i].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}

	for(i=0;i<totalnum.cresentpart;i++)
	{
		darea=ldcresentpartarea(starradius,i)*ellipse[cresentpart[i].ellipsen].fracdarkness;
#		if PRINTEACHSPOT
			ldspotreport[ellipse[cresentpart[i].ellipsen].spot]+=darea;
#		endif
		totarea-=darea;
	}

	if(totarea<0)
		if(totarea+ACCEPTABLEERROR>=0)
			totarea=0.0;
		else
			if(!errorflag)
			{
				//printf("\n(error 9 - area = %0.12lf <%2.6lf>)\n",totarea,totarea*1000000.0);	
				//fprintf(outerr,"\n(error 9 - area = %0.12lf <%2.6lf>)\n",totarea,totarea*1000000.0);
				//fflush(outerr);
				errorflag=9;
			}

	if(totarea>PI*starradius*starradius)
		if(!errorflag)
			errorflag=9;

	return totarea;
}
double fixldstararea(stardata *star,int ringi)
{
	double div,r0,r1,dr0,dr1,area0,area1;
	
	errorflag=0;
	if(ringi==0)
		r0=0;
	else
		r0=star->ringr[ringi-1];
	r1=star->ringr[ringi+1];
	dr0=star->ringr[ringi]-r0;
	dr1=r1-star->ringr[ringi];
	
	div=65536.0;
	r0=star->ringr[ringi]-dr0/div;
	r1=star->ringr[ringi]+dr1/div;
	area0=ldstararea(r0);
	area1=ldstararea(r1);
	while(errorflag&&div>=2.0)
	{
		errorflag=0;
		div/=2.0;
		r0=star->ringr[ringi]-dr0/div;
		r1=star->ringr[ringi]+dr1/div;
		area0=ldstararea(r0);
		area1=ldstararea(r1);
	}
	if(errorflag)
	{
		return 0.0;
	}

	return (area0+area1)/2.0;
}
double planetfindeccentricanomaly(double ma,double e)
{
	int i;
	double a,b,c;
	if(ma<PI)
	{
		a=0;
		c=PI;
	}
	else if(ma>PI)
	{
		a=PI;
		c=PIt2;
	}
	else
		return PI;
	for(i=0;i<16;i++)
	{
		b=(a+c)/2.0;
		if(b-e*sin(b)>ma)
			c=b;
		else
			a=b;
	}
	return (a+c)/2.0;
}
void planetposef(double t,planetdata planet[MAXPLANETS],double *x,double *y,int whichplanet)
{
	//planet's position in ellipse-of-orbit frame, x-axis is semimajor, origin at star, +x towards perastrion
	double ma,ea,cosea,costh;
	double theta,r;
	double xe,ye;
	
	ma=planet[whichplanet].omegaorbit*(t+planet[whichplanet].tmiddleoftransit);
	while(ma>=PIt2)
		ma-=PIt2;
	while(ma<0.0)
		ma+=PIt2;
	ea=planetfindeccentricanomaly(ma,planet[whichplanet].eccentricity);
	cosea=cos(ea);

	r=planet[whichplanet].orbitsemimajor*(1.0-planet[whichplanet].eccentricity*cosea);
	costh=(cosea-planet[whichplanet].eccentricity)/(1.0-planet[whichplanet].eccentricity*cosea);
	theta=acos(costh);

	xe=r*costh;
	ye=r*sin(theta);
	if(ma>PI) 
		ye=(-ye);

	(*x)=xe;
	(*y)=ye;

}
void setplanetpos(double t,planetdata planet[MAXPLANETS],int whichplanet)
{
	double xe,ye;

	planetposef(t,planet,&xe,&ye,whichplanet);
	planet[whichplanet].x=xe*planet[whichplanet].xhat[0]+ye*planet[whichplanet].yhat[0];
	planet[whichplanet].y=xe*planet[whichplanet].xhat[1]+ye*planet[whichplanet].yhat[1];
	planet[whichplanet].z=xe*planet[whichplanet].xhat[2]+ye*planet[whichplanet].yhat[2];
#	if XYZDETAILS
	 fprintf(xyzdetail,"%lf %lf %lf %lf\n",t,planet[whichplanet].x,planet[whichplanet].y,planet[whichplanet].z);
#	endif
}
void planetsetmiddleoftransit(planetdata planet[MAXPLANETS],int whichplanet)
{
	int i;
	double bt[2],bd[2],mt,md;
	double tmax,dt;

	planet[whichplanet].tmiddleoftransit=0.0;
	tmax=PIt2/planet[whichplanet].omegaorbit;
	dt=tmax/64.0;

	bd[0]=TOOBIG;
	bd[1]=TOOBIG;

	while(bd[1]==TOOBIG)
	{
		bd[0]=TOOBIG;
		bd[1]=TOOBIG;
		for(mt=0;mt<tmax;mt+=dt)
		{
			setplanetpos(mt,planet,whichplanet);
			if(planet[whichplanet].x>0)
			{
				md=planet[whichplanet].y*planet[whichplanet].y+planet[whichplanet].z*planet[whichplanet].z;
				if(md<bd[1])
				{
					if(md<bd[0])
					{
						bd[1]=bd[0];
						bt[1]=bt[0];
						bd[0]=md;
						bt[0]=mt;
					}
					else
					{
						bd[1]=md;
						bt[1]=mt;
					}
				}
			}
		}
		dt=dt/2.0;
	}


	if(bt[0]==0&&bt[1]>tmax-1.5*dt)
		bt[0]=tmax;
	if(bt[1]==0&&bt[0]>tmax-1.5*dt)
		bt[1]=tmax;

	for(i=0;i<64;i++)
	{
		mt=(bt[0]+bt[1])/2.0;
		setplanetpos(mt,planet,whichplanet);
		md=planet[whichplanet].y*planet[whichplanet].y+planet[whichplanet].z*planet[whichplanet].z;
		if(bd[0]<bd[1])
		{
			bd[1]=md;
			bt[1]=mt;
		}
		else
		{
			bd[0]=md;
			bt[0]=mt;
		}
	}
	planet[whichplanet].tmiddleoftransit=(bt[0]+bt[1])/2.0;
}

double lightness(double t,stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS])
{
	int i,j,realnumplanets;
	double totlight,morelight,torig,thelightness,noplanetlightness;

	torig=t/86400.0+planet->lct0;
#	if ANYPRINTVIS
		if(PRINTVIS)
			if(PRINTVIS==1)
				fprintf(outv,"%.9lf %10.6lf %10.6lf %10.6lf %10.6lf 9\n",torig,0.0,0.0,1.0,1.0);
			else
			{
				fprintf(outv,"c\n0 0 %lf 9\n",star->ringr[NLDRINGS-1]); 
			}
#	endif
	zerototalnums();

	star->phi=star->omega*t;

	for(i=0;i<numplanets;i++)
	{
		setplanetpos(t,planet,i);
		if(planet[i].x>0)	//planet is in front of star
		{
			if(sqrt(planet[i].y*planet[i].y+planet[i].z*planet[i].z)-planet[i].r<star->r)
			{
				createcircle(planet[i],star->r,star->rsq);
#				if	ANYPRINTVIS
					if(PRINTVIS)
						if(PRINTVIS==1)
							fprintf(outv,"%.9lf %10.6lf %10.6lf %10.6lf %10.6lf 3\n",torig,planet[i].y,planet[i].z,planet[i].r,planet[i].r);
						else
							fprintf(outv,"c\n%g %g %g 3\n",planet[i].y,planet[i].z,planet[i].r);
							
#				endif
			}
#			if ANYPRINTVIS
				else if(PRINTVIS)
					if(PRINTVIS==1)
						fprintf(outv,"%.9lf %10.6lf %10.6lf %10.6lf %10.6lf 10\n",torig,planet[i].y,planet[i].z,planet[i].r,planet[i].r);
					else
						fprintf(outv,"c\n%g %g %g 10\n",planet[i].y,planet[i].z,planet[i].r);
						
#			endif
		}
	}

	updatespots(spot,star->r,star->phi);
	for(i=0;i<numspots;i++)
		if(spot[i].psi-spot[i].alpha<PIo2)	//spot visable
			createellipsecresent(i,spot[i],star->r);
#	if ANYPRINTVIS
		if(PRINTVIS)
			if(PRINTVIS==1)
			{
				for(i=0;i<numspots;i++)
				{
					if(spot[i].psi>PIo2+spot[i].alpha)
						fprintf(outv,"%.9lf %10.6lf %10.6lf %10.6lf %10.6lf 10\n",torig,spot[i].y,spot[i].z,spot[i].r,spot[i].r*spot[i].cospsi);
					else
						fprintf(outv,"%.9lf %10.6lf %10.6lf %10.6lf %10.6lf 1\n",torig,spot[i].y,spot[i].z,spot[i].r,spot[i].r*spot[i].cospsi);

				}
			}
			else
				for(i=0;i<numspots;i++)
					if(spot[i].psi-spot[i].alpha>PIo2)
						fprintf(outv,"e\n%g %g %g %g %i\n",spot[i].y,spot[i].z,spot[i].r,spot[i].r*spot[i].cospsi,pvisc);
					else
						fprintf(outv,"e\n%g %g %g %g %i\n",spot[i].y,spot[i].z,spot[i].r,spot[i].r*spot[i].cospsi,pvasc);
#	endif					
	checkoverlap(star->r);	//create disjoint shapes that block light

	totlight=realstararea(star->r)*star->dintensity[NLDRINGS-1];
	if(PRINTEACHSPOT) 
		for(i=0;i<numspots;i++)
			spotreport[i]=ldspotreport[i]*star->dintensity[NLDRINGS-1];

	for(i=NLDRINGS-2;i>=0;i--)
	{
		debugringi=i;
		morelight=ldstararea(star->ringr[i])*star->dintensity[i];
		if(errorflag>100)
			morelight=fixldstararea(star,i)*star->dintensity[i];
/*		if(errorflag)
		{
			printf("error %i\n",errorflag);	
			fprintf(outerr,"error trial %i lcni %i ringi %i errornum %i\n",debugtrial,debuglcni,debugringi,errorflag);
			fflush(outerr);
			errorflag=0;
		}*/
		if(PRINTEACHSPOT)
			for(j=0;j<numspots;j++)
				spotreport[j]+=ldspotreport[j]*star->dintensity[i];
		totlight+=morelight;
	}

	if(PRINTEACHSPOT)
		for(i=0;i<numspots;i++)
			spotreport[i]*=star->brightnessfactor*star->area/star->maxlight;
	if(!FLATTENMODEL||numplanets==0)
		thelightness=star->brightnessfactor*star->area*totlight/star->maxlight;	//normalized to PI*Rstar^2 * brightnessfactor
	else
	{
		thelightness=star->brightnessfactor*star->area*totlight/star->maxlight;
		realnumplanets=numplanets;
		numplanets=0;
		noplanetlightness=lightness(t,star,planet,spot);
		thelightness=thelightness-noplanetlightness+PI;  //flattened to remove spot effects except planet-spot crossings
		numplanets=realnumplanets;
	}
	return thelightness;
}
#if PRECALCPLANET
double lightnessindexed(double t,stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int timeindex)
{
	int i,j,realnumplanets;
	double totlight,morelight,torig,thelightness,noplanetlightness;

#	if ANYPRINTVIS
		if(PRINTVIS)
			return lightness(t,star,planet,spot);	//doesn't work with printvis
#	endif

	torig=t/86400.0+planet->lct0;
	zerototalnums();

	star->phi=star->omega*t;

	totalnum.circle=totalnumcircle[timeindex];
	for(i=0;i<totalnum.circle;i++)
	{
		circle[i].centery=indcircle[i][timeindex].centery; circle[i].centerz=indcircle[i][timeindex].centerz;
		circle[i].radius=indcircle[i][timeindex].radius; circle[i].area=indcircle[i][timeindex].area; circle[i].rsq=indcircle[i][timeindex].rsq;
		circle[i].mindcen=indcircle[i][timeindex].mindcen;circle[i].maxdcen=indcircle[i][timeindex].maxdcen;
		circle[i].numstarint=indcircle[i][timeindex].numstarint;
		if(circle[i].numstarint>0)
		{
			circle[i].starintdang=indcircle[i][timeindex].starintdang;
			circle[i].starintbeta[0]=indcircle[i][timeindex].starintbeta[0];circle[i].starintbeta[1]=indcircle[i][timeindex].starintbeta[1];
		}
	}

	updatespots(spot,star->r,star->phi);
	for(i=0;i<numspots;i++)
		if(spot[i].psi-spot[i].alpha<PIo2)	//spot visable
			createellipsecresent(i,spot[i],star->r);
	checkoverlap(star->r);	//create disjoint shapes that block light

	totlight=realstararea(star->r)*star->dintensity[NLDRINGS-1];
	if(PRINTEACHSPOT) 
		for(i=0;i<numspots;i++)
			spotreport[i]=ldspotreport[i]*star->dintensity[NLDRINGS-1];

	for(i=NLDRINGS-2;i>=0;i--)
	{
		debugringi=i;
		morelight=ldstararea(star->ringr[i])*star->dintensity[i];
		if(errorflag>100)
			morelight=fixldstararea(star,i)*star->dintensity[i];
		if(PRINTEACHSPOT)
			for(j=0;j<numspots;j++)
				spotreport[j]+=ldspotreport[j]*star->dintensity[i];
		totlight+=morelight;
	}

	if(PRINTEACHSPOT)
		for(i=0;i<numspots;i++)
			spotreport[i]*=star->brightnessfactor*star->area/star->maxlight;
	if(!FLATTENMODEL||numplanets==0)
		thelightness=star->brightnessfactor*star->area*totlight/star->maxlight;	//normalized to PI*Rstar^2 * brightnessfactor
	else
	{
		thelightness=star->brightnessfactor*star->area*totlight/star->maxlight;
		realnumplanets=numplanets;
		numplanets=0;
		noplanetlightness=lightness(t,star,planet,spot);
		thelightness=thelightness-noplanetlightness+PI;  //flattened to remove spot effects except planet-spot crossings
		numplanets=realnumplanets;
	}
	return thelightness;
}
#endif
int inifilereaderld(char datastr[4096],int *a,int e,double lda[5])
{
	int i,b[4],c;

	if(datastr[*a]>57||datastr[*a]<45||datastr[*a]==47)
	{
		printf("line doesn't start with a number\n");
		return -1;
	}

	b[0]=*a;
	for(i=0;i<3;i++)
	{
		while(datastr[b[i]]!=32&&datastr[b[i]]!=9)
		{
			b[i]++;
			if(b[i]>=e)
				return -1;
		}
		datastr[b[i]]=0;
		b[i+1]=b[i]+1;
	}
	while(datastr[b[3]]!=32&&datastr[b[3]]!=9&&datastr[b[3]]!=10&&datastr[b[3]]!=0)
	{
		b[3]++;
		if(b[3]>=e)
			return -1;
	}
	c=b[3];
	if(datastr[c]!=0)
		while(datastr[c]!=0)
		{
			c++;
			if(c>=e)
				return -1;
		}
	while(datastr[c]==0)
	{
		c++;
		if(c>=e)
			return -1;
	}
	datastr[b[3]]=0;
	sscanf(datastr+*a,"%lf",&lda[1]);
	sscanf(datastr+b[0]+1,"%lf",&lda[2]);
	sscanf(datastr+b[1]+1,"%lf",&lda[3]);
	sscanf(datastr+b[2]+1,"%lf",&lda[4]);

#	if !QUIET
		printf("limb darkening: %lf %lf %lf %lf\n",lda[1],lda[2],lda[3],lda[4]);
#	endif

	*a=c;
	return 0;
}
double ldintensity(double r,double lda[5])
{	//intensity at radius r (limb darkening)
	double mu;
	mu=sqrt(1-r*r);
	return 1-lda[1]*(1-sqrt(mu))-lda[2]*(1-mu)-lda[3]*(1-pow(mu,1.5))-lda[4]*(1-mu*mu);
}
double ldradiusfromintensity(double intensity,double lda[5])
{
	int count=0;
	double minr=0.0,maxr=1.0,iatmin,iatmax;
	double midr,iatmid,di;
	double rfitollerance=0.0001;

	iatmin=ldintensity(minr,lda);
	iatmax=ldintensity(maxr,lda);
	di=100;

	if(iatmin<intensity||iatmax>intensity)
		return -1.0;

	while(di>rfitollerance)
	{
		count++;
		midr=(minr+maxr)/2.0;
		iatmid=ldintensity(midr,lda);
		if(iatmid>intensity)
		{
			minr=midr;
			di=iatmid-intensity;
		}
		else
		{
			maxr=midr;
			di=intensity-iatmid;
		}
		if(count>100000000)
			di=0;
	}

	return midr;
}
double ldmaxintensity(double lda[5])
{
	int i;
	double rmin,rmid,rmax,imin,imid,imax;

	imin=0.0;
	imax=0.0;
	rmax=0.0;
	for(rmid=0.0;rmid<=1.0;rmid+=0.001)
	{
		imid=ldintensity(rmid,lda);
		if(imid>imin)
			if(imid>imax)
			{
				rmin=rmax;
				imin=imax;
				rmax=rmid;
				imax=imid;
			}
			else
			{
				rmin=rmid;
				imin=imid;
			}
	}

	if(rmin>rmax)
	{
		rmid=rmin;
		rmin=rmax;
		rmax=rmid;
	}
	imin=ldintensity(rmin,lda);
	imax=ldintensity(rmax,lda);
	for(i=0;i<100;i++)
	{
		rmid=(rmin+rmax)/2.0;
		if(imin>imax)
		{
			rmax=rmid;
			imax=ldintensity(rmax,lda);
		}
		else if(imin<imax)
		{
			rmin=rmid;
			imin=ldintensity(rmin,lda);
		}
		else
		{
			i=100;
		}
	}

	if(imax>imin)
		return imax;
	else
		return imin;
}
int initldrings(stardata *star,double lda[5])
{
	int i;
	double delta,imax;
	double totintensity[MAXNLDRINGS];

	if(NLDRINGS==1)
	{
		star->ringr[0]=star->r;
		star->dintensity[0]=1.0;
		return 1;
	}

#	if USEOLDLDRING
		for(i=0;i<NLDRINGS;i++)
			star->ringr[i]=(i+1)*star->r/NLDRINGS; //outer radius of ith ring
		star->dintensity[0]=ldintensity(star->ringr[0]/2.0,lda)-ldintensity((star->ringr[0]+star->ringr[1])/2.0,lda);
		for(i=1;i<NLDRINGS-1;i++)	//amount intensity of ring i is greater than intensity of ring i+1
			star->dintensity[i]=ldintensity((star->ringr[i-1]+star->ringr[i])/2.0,lda)-ldintensity((star->ringr[i]+star->ringr[i+1])/2.0,lda);
		star->dintensity[NLDRINGS-1]=ldintensity((star->ringr[NLDRINGS-2]+star->ringr[NLDRINGS-1])/2.0,lda);
		return 1;
#	else
		imax=ldmaxintensity(lda);
		delta=(imax-ldintensity(1.0,lda))/NLDRINGS;

		for(i=0;i<NLDRINGS-1;i++)
			star->ringr[i]=ldradiusfromintensity(imax-((double)i+1.0)*delta,lda);
		star->ringr[NLDRINGS-1]=1.0;
		totintensity[0]=ldintensity(star->ringr[0]/2.0,lda);
		star->dintensity[0]=ldintensity(star->ringr[0]/2.0,lda)-ldintensity((star->ringr[0]+star->ringr[1])/2.0,lda);
		for(i=1;i<NLDRINGS-1;i++)	//amount intensity of ring i is greater than intensity of ring i+1
		{
			totintensity[i]=ldintensity((star->ringr[i-1]+star->ringr[i])/2.0,lda);
			star->dintensity[i]=ldintensity((star->ringr[i-1]+star->ringr[i])/2.0,lda)-ldintensity((star->ringr[i]+star->ringr[i+1])/2.0,lda);
		}
		totintensity[NLDRINGS-1]=ldintensity((star->ringr[NLDRINGS-2]+star->ringr[NLDRINGS-1])/2.0,lda);
		star->dintensity[NLDRINGS-1]=ldintensity((star->ringr[NLDRINGS-2]+star->ringr[NLDRINGS-1])/2.0,lda);
		return 1;
#	endif
}
int prepfileread(char fn[128],char *datastr,int maxlength)
{	//read file into string, return length of string
	char ch;
	int i;
	FILE *in;

	in=fopen(fn,"r");
	if(in==NULL)
		return 0;
	i=0;
	while(!feof(in))
	{
		ch=fgetc(in);
		if(ch==13)
			ch=10;	//CR  to LF
		if(ch==9||ch==44)
			ch=32;	//Tab and , to Space
		if(ch==(-1))
			ch=0;
		if(i==0&&ch!=10&&ch!=32)
		{
			datastr[i]=ch;
			i++;
		}
		else if(i!=0)
		{
			if(ch!=10&&ch!=32)
			{
				datastr[i]=ch;
				i++;
			}
			else
			{
				if(ch==32)
				{
					if(datastr[i-1]!=32&&datastr[i-1]!=0)
					{
						datastr[i]=32;
						i++;
					}
				}
				else if(ch==10)
				{
					if(datastr[i-1]==32)
						datastr[i-1]=0;
					else if(datastr[i-1]!=0)
					{
						datastr[i]=0;
						i++;
					}
				}
			}
		}
		if(i>=maxlength)
			return (-1);
	}
	if(datastr[i-1]!=0)
	{
		datastr[i]=0;
		i++;
	}
	return i;
	fclose(in);
}
int filereadd(int n,double *x,char *datastr,int *current,int end)
{	//reads n doubles from fn and puts them in x[] return 0 if succesful, -1 if not
	int i,a,b;
	double y;

	a=*current;
	for(i=0;i<n;i++)
	{
		if(datastr[a]>57||datastr[a]<45||datastr[a]==47)
		{
			printf("input is not a number\n");
			return -1;
		}
		b=a;
		while((datastr[b]>=48&&datastr[b]<=57)||datastr[b]==46||datastr[b]==45)
		{
			b++;
			if(b>=end)
			{
				printf("not enough numbers in file\n");
				return -1;
			}
		}
		if(b==end-1)
			if(i==n-1)
			{
				sscanf(datastr+a,"%lf",&y);
				x[i]=y;
#				if !QUIET
					printf("(%lf)\n",x[i]);
#				endif
				return 0;
			}
			else
				return -1;

		if(datastr[b]!=0)
		{
			datastr[b]=0;
			sscanf(datastr+a,"%lf",&y);
			x[i]=y;
			b++;
			while(datastr[b]==32&&b<end)
				b++;
#			if !QUIET
				printf("(%lf) \"%s\"\n",x[i],datastr+b);
#			endif
		}
		else
		{
			sscanf(datastr+a,"%lf",&y);
			x[i]=y;
#			if !QUIET
				printf("(%lf)\n",x[i]);
#			endif
		}

		while(datastr[b]!=0)
		{
			b++;
			if(b>=end)
				return -1;
		}
		while(datastr[b]==0)
		{
			b++;
			if(b>=end)
				if(i<n-1)
					return -1;
				else
				{
					*current=b;
					return 0;
				}
		}

		a=b;
		if(i==n-1)
		{
			*current=a;
			return 0;
		}
	}
	return -1;
}
int filereads(char *rstr,char *datastr,int *current,int end)
{	//reads one line from datastr and puts it into rstr. return 0 if succesful, -1 if not
	int a;

	a=*current;
	sscanf(datastr+a,"%s",rstr);
#	if !QUIET
		printf("\'%s\'\n",rstr);
#	endif
	a++;
	while(datastr[a]!=0)
	{
		a++;
		if(a>=end)
			return -1;
	}
	a++;
	*current=a;
	return 0;
}

int initializestarplanet(stardata *star,planetdata planet[MAXPLANETS],char filename[64],double *lcstarttime,double *lcfinishtime,double *lcmaxlight,double *ascale,int *mcmcnpop,long int *mcmcmaxstepsortime,int *partitionpop,int *partitionsteps,double *readparam,int *randomseed,char seedfilename[64])
{
	int i,a,b,e;
	char str[128],datastr[4096];
	double x[9],stardensity,ppdays;
	double lda[5];		//limb darkening coeffecients (really just 4)
	
	e=prepfileread(filename,datastr,4096);
	if(e==0)
		return -1;
	else if(e<0)
		return -2;

	a=0;
	if(filereads(str,datastr,&a,e)<0)
		return -3;
	if(str[0]!='#')	//planet properties
		return -3;

	if(filereadd(1,x,datastr,&a,e)<0)
		return -4;
	numplanets=(int)x[0];
#	if !QUIET
		printf("number of planets (%i)\n",numplanets);
#	endif
	if(numplanets>MAXPLANETS)
	{
		printf("too many planets\n");
		return -5;
	}

	planet[0].lct0=0.0;	//in case numplanets==0
	for(i=0;i<numplanets;i++)
	{
		if(filereadd(9,x,datastr,&a,e)<0)
			return -6;
		planet[i].lct0=x[0];
		planet[i].omegaorbit=PIt2/(x[1]*86400.0);
		ppdays=x[1];
		planet[i].rsq=x[2]; //needs to be multipled by star->rsq
		x[5]=x[5]*PI/180.0;
		x[6]=x[6]*PI/180.0;
		planet[i].thetaorbit=acos(cos(x[6])*sin(x[5]));
		planet[i].phiorbit=atan2((-1.0)*sin(x[6])*sin(x[5]),cos(x[5]));
		if(planet[i].thetaorbit<0)
			planet[i].thetaorbit*=(-1.0);
		if(planet[i].phiorbit<0)
			planet[i].phiorbit+=PIt2;
		if(x[7]==0&&x[8]==0)
		{
			planet[i].orbitangleomega=PIo2;
			planet[i].eccentricity=0;
		}
		else if(x[7]==0)
		{
			planet[i].orbitangleomega=PIo2;
			planet[i].eccentricity=x[8];
		}
		else
		{
			planet[i].orbitangleomega=atan2(x[8],x[7]);
			planet[i].eccentricity=x[7]/cos(planet[i].orbitangleomega);
		}
//		planet[i].orbitangleomega+=PI;  //not PIo2;	//because 90 degrees is star at far focus, not near focus 
//		if(planet[i].orbitangleomega>=PIt2)
//			planet[i].orbitangleomega-=PIt2;	perhaps we thought better of this?

#		if !QUIET
			printf("planet %i\n orbit theta= %0.9lf (%lf degrees)\n orbit phi= %0.9lf (%lf degrees)\n",i,planet[i].thetaorbit,planet[i].thetaorbit*180.0/PI,planet[i].phiorbit,planet[i].phiorbit*180/PI);
			printf(" eccentricity= %lf\n orbit angle omega= %lf (%lf degrees)\n",planet[i].eccentricity,planet[i].orbitangleomega,planet[i].orbitangleomega*180.0/PI);
#		endif
/*		planet[i].thetaorbit=90.0-x[5];
		if(planet[i].thetaorbit<0)
			planet[i].thetaorbit*=(-1);
		planet[i].thetaorbit*=PI/180.0;
		planet[i].phiorbit=0.0;  */
	}
	
	if(filereads(str,datastr,&a,e)<0)
		return -7;
	if(str[0]!='#')		//star properties
		return -7;

	if(filereadd(5,x,datastr,&a,e)<0)
		return -8;
	stardensity=x[0];
	if(x[1]!=0.0)
		star->omega=PIt2/(x[1]*86400.0);
	else
		star->omega=0.0;
	star->theta=x[4]*PI/180.0;

	if(inifilereaderld(datastr,&a,e,lda)<0)
		return -9;

	if(filereadd(1,x,datastr,&a,e)<0)
		return -9;
	NLDRINGS=(int)x[0];
	if(NLDRINGS>MAXNLDRINGS)
	{
		printf("too many rings\n");
		return -9;
	}

	if(filereads(str,datastr,&a,e)<0)
		return -10;
	if(str[0]!='#')		//spot properties
		return -10;

	if(filereadd(2,x,datastr,&a,e)<0)
		return -11;
	numspots=(int)x[0];
#	if !QUIET
		printf("number of spots (%i)\n",numspots);
#	endif
	if(numspots>MAXSPOTS)
	{
		printf("too many spots\n");
		return -12;
	}
	spotfracdark=1.0-x[1];
	
	if(filereads(str,datastr,&a,e)<0)
		return -13;
	if(str[0]!='#')		//light curve
		return -13;

	if(filereads(filename,datastr,&a,e)<0)
		return -14;
#	if !QUIET
		printf("filename for lightcurve data (%s)\n",filename);
#	endif
	if(filereadd(4,x,datastr,&a,e)<0)
		return -15;
	*lcstarttime=x[0];
	*lcfinishtime=x[0]+x[1];
	*lcmaxlight=x[2];
	if((*lcmaxlight)==0)
	{
		USEDOWNFROMMAX=1;
#		if !QUIET
			printf("using down from max (%i)\n",DOWNFROMMAX);
#		endif
	}
	else
		USEDOWNFROMMAX=0;
	FLATTENMODEL=(int)x[3];
#	if !QUIET
		if(FLATTENMODEL)
			printf("flattening model light curve\n");
		else
			printf("not flattening model light curve\n");
#	endif
	star->r=1.0;
	if(!initldrings(star,lda))
		return -31;
	star->maxlight=0.0;
	for(i=NLDRINGS-1;i>=0;i--)
		star->maxlight+=PI*star->ringr[i]*star->ringr[i]*star->dintensity[i];	
	star->rsq=star->r*star->r;
	star->area=PI*star->rsq;
	star->brightnessfactor=1.0;

	for(i=0;i<numplanets;i++)
	{
		double cw90,sw90,st,ct,sf,cf,p;
		planet[i].rsq*=star->rsq;
		planet[i].r=sqrt(planet[i].rsq);
		planet[i].area=PI*planet[i].rsq;
		planet[i].orbitsemimajor=pow(ppdays/365.24218967,2.0/3.0)*pow(stardensity,1.0/3.0)*214.939469384;
#		if !QUIET
			printf("planet %i orbit semimajor axis= %lf\n",i,planet[i].orbitsemimajor);
#		endif
		cw90=cos(planet[i].orbitangleomega-PIo2);
		sw90=sin(planet[i].orbitangleomega-PIo2);
		ct=cos(planet[i].thetaorbit);
		st=sin(planet[i].thetaorbit);
		cf=cos(planet[i].phiorbit);
		sf=sin(planet[i].phiorbit);
		p=sqrt(st*st*sf*sf+ct*ct);
		planet[i].xhat[0]=cw90*p;
		planet[i].xhat[1]=(sw90*ct-cw90*st*st*sf*cf)/p;
		planet[i].xhat[2]=(-cw90*st*ct*cf-sw90*st*sf)/p;
		planet[i].yhat[0]=(-sw90*p);
		planet[i].yhat[1]=(sw90*st*st*sf*cf+cw90*ct)/p;
		planet[i].yhat[2]=(sw90*st*ct*cf-cw90*st*sf)/p;
		planetsetmiddleoftransit(planet,i);
	}

	if(filereads(str,datastr,&a,e)<0)
		return -16;
	if(str[0]!='#')		//action 
		return -16;

	if(filereads(str,datastr,&a,e)<0)
		return -17;

	if(str[0]!='L'&&str[0]!='l'&&str[0]!='H'&&str[0]!='h'&&str[0]!='S'&&str[0]!='s'&&str[0]!='D'&&str[0]!='d'&&str[0]!='T'&&str[0]!='t'&&str[0]!='x'&&str[0]!='X'&&str[0]!='f'&&str[0]!='F'&&str[0]!='p'&&str[0]!='P'&&str[0]!='u'&&str[0]!='U'&&str[0]!='I'&&str[0]!='i'&&str[0]!='a'&&str[0]!='A')
		i=0;  //unseeded mcmc
	else if(str[0]=='l')
		i=1;	//generate light curve and vis file from parameters
	else if(str[0]=='H'||str[0]=='h')
		i=2;	//metropolis-hastings
	else if(str[0]=='s')
		i=3;	//single seeded mcmc from parameters
	else if(str[0]=='D'||str[0]=='d')
		i=4;	//debug
	else if(str[0]=='T'||str[0]=='t')
		i=5;	//totally seeded mcmc
	else if(str[0]=='L')
		i=6;	//generate light curve and vis file from parameter file
	else if(str[0]=='S')
		i=7;	//single seeded mcmc from parameter file
	else if(str[0]=='f'||str[0]=='F')
	{
		FIXTHETAS=1;
		if(str[1]=='s')
			i=3;
		else if(str[1]=='S')
			i=7;
		else if(str[1]=='T'||str[1]=='t')
			i=5;
		else
			i=8;	//unseeded mcmc with fixed theta
	}
	else if(str[0]=='p'||str[0]=='P')
		i=9;	//plot chi squared over variation of one spot
	else if(str[0]=='u'||str[0]=='U')
		i=10;	//partially seeded mcmc
	else if(str[0]=='i'||str[0]=='I')
		i=11;
	else if(str[0]=='a')
		i=12;
	else if(str[0]=='A')
		i=13;
	else if(str[0]=='x')
		i=100;	//x
	else
		i=101;	//X

	numseeded=0;  //unless set below
	if(i==0||i==3||i==5||i==7||i==8||i==10)
	{
		if(filereadd(5,x,datastr,&a,e)<0)
			return -18;
		*randomseed=(int)x[0];
		*ascale=x[1];
		*mcmcnpop=(int)x[2];
		*mcmcmaxstepsortime=(long int)x[3];
		CALCBRIGHTNESSFACTOR=(char)x[4];
		if(i==3||i==7)
		{
			if(filereadd(2,x,datastr,&a,e)<0)
				return -19;
			sigmaradius=x[0];
			sigmaangle=x[1];
		}
		if(i==5||i==7)
		{
			if(filereads(seedfilename,datastr,&a,e)<0)
				return -20;
			if(i==7)
			{
				FILE *in;
				printf("parameter file: %s\n",seedfilename);
				in=fopen(seedfilename,"r");
				if(in==NULL)
					return -21;
				for(b=0;b<numspots*3+1;b++)
					fscanf(in,"%lf\n",&readparam[b]);
				fclose(in);
			}
		}
		else if(i==3)
		{
			if(filereadd(numspots*3+1,readparam,datastr,&a,e)<0)
				return -22;
		}
		else if(i==8)
		{
			if(filereadd(numspots,readparam,datastr,&a,e)<0)
				return -29;

		}
		else if(i==10)
		{
			if(filereadd(2,x,datastr,&a,e)<0)
				return -30;
			if(x[0]>=numspots)
				return -31;
			numseeded=(int)x[0];	//number of seeded spots
			FIXSEEDEDONLYPHI=(int)x[1];
			if(FIXSEEDEDONLYPHI)         //only phi or fix none
			{   
				if(filereadd(2,x,datastr,&a,e)<0)
					return -32;
				sigmaradius=x[0];
				sigmaangle=x[1];
			}
			if(filereadd(numseeded*3+1,readparam,datastr,&a,e)<0)
				return -33;
		}
	}
	else if(i==1||i==12)
	{	
		if(filereadd(numspots*3+1,readparam,datastr,&a,e)<0)
				return -23;
	}
	else if(i==6||i==13)
	{
		if(filereads(seedfilename,datastr,&a,e)<0)
			return -24;
		FILE *in;
		printf("parameter file: %s\n",seedfilename);
		in=fopen(seedfilename,"r");
		for(b=0;b<numspots*3+1;b++)
			fscanf(in,"%lf\n",&readparam[b]);
		fclose(in);
		b=0;
		while(seedfilename[b]!=0&&b<64)
			b++;
		if(seedfilename[b]!=0)
			return -25;
		b-=13;
		if(seedfilename[b]=='p'&&seedfilename[b+1]=='a'&&seedfilename[b+2]=='r'&&seedfilename[b+3]=='a'&&seedfilename[b+4]=='m'&&seedfilename[b+5]=='b'&&seedfilename[b+6]=='e'&&seedfilename[b+7]=='s'&&seedfilename[b+8]=='t'&&seedfilename[b+9]=='.'&&seedfilename[b+10]=='t'&&seedfilename[b+11]=='x'&&seedfilename[b+12]=='t')
		{
			seedfilename[b]='t'; seedfilename[b+1]='v'; seedfilename[b+2]='i'; seedfilename[b+3]='s'; seedfilename[b+4]='.'; seedfilename[b+5]='t'; seedfilename[b+6]='x'; seedfilename[b+7]='t'; seedfilename[b+8]=0;
		}
		else
		{
			b+=13;
			seedfilename[b]='_';
			b++;
			seedfilename[b]='t'; seedfilename[b+1]='v'; seedfilename[b+2]='i'; seedfilename[b+3]='s'; seedfilename[b+4]='.'; seedfilename[b+5]='t'; seedfilename[b+6]='x'; seedfilename[b+7]='t'; seedfilename[b+8]=0;
		}
	}
	else if(i==2)
	{
		if(filereadd(3,x,datastr,&a,e)<0)
			return -26;
		*mcmcmaxstepsortime=(long int)x[0];
		sigmaradius=x[1];
		sigmaangle=x[2];
		if(filereadd(numspots*3+1,readparam,datastr,&a,e)<0)
			return -27;
	}
	else if(i==9)
	{
		if(filereadd(1,x,datastr,&a,e)<0)
			return -28;
		*randomseed=(int)x[0]; //using randomseed to hold number of spot to vary
		if(filereadd(numspots*3+1,readparam,datastr,&a,e)<0)
			return -29;
	}
	else if(i==11)
	{
		if(filereads(seedfilename,datastr,&a,e)<0)
			return -30;
		printf("time file: %s\n",seedfilename);

	}
	return i;
}
int lcreadline(FILE *in,double x[3])
{	//returns -1 if hits eof
	char datastr[128],ch;
	int i;
	double a,b,c;

	ch=10;
	while(ch==10||ch==13)
	{
		ch=fgetc(in);
		if(feof(in))
			return -1;
	}
	while(ch=='#')
	{
		while(ch!=10&&ch!=13)
		{
			ch=fgetc(in);
			if(feof(in))
				return -1;
		}
		while(ch==10||ch==13)
		{
			ch=fgetc(in);
			if(feof(in))
				return -1;
		}
	}
	datastr[0]=ch;
	i=0;
	while(ch!=10&&ch!=13)
	{
		ch=fgetc(in);
		i++;
		datastr[i]=ch;
		if(feof(in))
			return -1;
	}
	datastr[i]=0;
	sscanf(datastr,"%lf %lf %lf",&a,&b,&c);
	x[0]=a;
	x[1]=b;
	x[2]=c;
	return 0;
}
void preinitializelcdata(char filename[64],double lcstarttime,double lcfinishtime,double lcmaxlight,int *lcn,double *lclightnorm)
{
	char go,inrange;
	int i,j,n;
	double x[3],maxl[DOWNFROMMAX];
	FILE *in;

	in=fopen(filename,"r");
	if(USEDOWNFROMMAX)
		for(i=0;i<DOWNFROMMAX;i++)
			maxl[i]=0;

	n=0;
	inrange=0;
	go=lcreadline(in,x);
	while(go>=0)
	{
		if(x[0]>=lcstarttime&&x[0]<=lcfinishtime)
			n++;
		if(USEDOWNFROMMAX)
			if(x[1]>maxl[DOWNFROMMAX-1])
				for(i=0;i<DOWNFROMMAX;i++)
					if(x[1]>maxl[i])
					{
						for(j=DOWNFROMMAX-1;j>i;j--)
							maxl[j]=maxl[j-1];
						maxl[i]=x[1];
						i=DOWNFROMMAX;
					}

		go=lcreadline(in,x);
	}

	(*lcn)=n;
	if(USEDOWNFROMMAX)
		(*lclightnorm)=PI/maxl[DOWNFROMMAX-1];
	else
		(*lclightnorm)=PI/lcmaxlight;
	fclose(in);
}
int initializelcdata(char filename[64],double lcstarttime,double lcfinishtime,int lcn,double lclightnorm,double lctime[],double lclight[],double lcuncertainty[],planetdata *planet)
{
	int i,go;
	double x[3];
	FILE *in;

	in=fopen(filename,"r");

	go=lcreadline(in,x);
	if(go<0)
		return -1;
	while(go>=0)
	{
		if(x[0]>=lcstarttime)
			go=(-1);
		else
		{
			go=lcreadline(in,x);
			if(go<0)
				return -2;
		}
	}

	for(i=0;i<lcn;i++)
	{
		lctime[i]=(x[0]-planet[0].lct0)*86400.0;	// if numplanets=0, planet[0].lct0=0
		lclight[i]=x[1]*lclightnorm;
		lcuncertainty[i]=x[2]*lclightnorm;
		if(i<lcn-1)
		{
			go=lcreadline(in,x);
			if(go<0)
				return -3;
		}
	}

	fclose(in);
	return 0;
	
}
#if PRECALCPLANET
void initializeprecalcplanet(stardata *star,planetdata planet[MAXPLANETS],int lcn,double lctime[])
{
	int i,timeindex;
	double torig;

	for(timeindex=0;timeindex<lcn;timeindex++)
	{
		torig=lctime[timeindex]/86400.0+planet->lct0;
		zerototalnums();

		for(i=0;i<numplanets;i++)
		{
			setplanetpos(lctime[timeindex],planet,i);
			if(planet[i].x>0)
			{
				if(sqrt(planet[i].y*planet[i].y+planet[i].z*planet[i].z)-planet[i].r<star->r)
				{
					createcircle(planet[i],star->r,star->rsq);
				}
			}
		}
		totalnumcircle[timeindex]=totalnum.circle;
		for(i=0;i<totalnumcircle[timeindex];i++)
		{
			indcircle[i][timeindex].centery=circle[i].centery; indcircle[i][timeindex].centerz=circle[i].centerz;
			indcircle[i][timeindex].radius=circle[i].radius; indcircle[i][timeindex].area=circle[i].area; indcircle[i][timeindex].rsq=circle[i].rsq;
			indcircle[i][timeindex].mindcen=circle[i].mindcen;indcircle[i][timeindex].maxdcen=circle[i].maxdcen;
			indcircle[i][timeindex].numstarint=circle[i].numstarint;
			if(indcircle[i][timeindex].numstarint>0)
			{
				indcircle[i][timeindex].starintdang=circle[i].starintdang;
				indcircle[i][timeindex].starintbeta[0]=circle[i].starintbeta[0]; indcircle[i][timeindex].starintbeta[1]=circle[i].starintbeta[1];
			}
		}
	}

}
#endif
double findchisq(stardata *star,planetdata *planet,spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[])
{
	int i;
	double dif,chisq;

	chisq=0;
	for(i=0;i<lcn;i++)
	{
		errorflag=0;
		debuglcni=i;
#		if PRECALCPLANET
			dif=lightnessindexed(lctime[i],star,planet,spot,i)-lclight[i];
#		else
			dif=lightness(lctime[i],star,planet,spot)-lclight[i];
#		endif
		if(errorflag)
		{
			fprintf(outerr,"error %i\nat datum %i\n (location 0)\n",errorflag,i);
			return -1;
		}
		chisq+=dif*dif/(lcuncertainty[i]*lcuncertainty[i]);
	}
	return chisq;
}
double normaldist(double x0,double sigma)
{
	double u,v,s;
	u=2.0*RND-1.0;
	v=2.0*RND-1.0;
	s=u*u+v*v;
	while(s>=u||s>=v)
	{
		u=2.0*RND-1.0;
		v=2.0*RND-1.0;
		s=u*u+v*v;
	}
	return x0+sigma*u*sqrt(-2.0*log(s)/s);
}
void setrandomparam(double p[],stardata *star)
{
	int i,j,k;
	double posunit[MAXSPOTS][3],alpha[MAXSPOTS],rspot,theta,phi,ang; //just to avoid overlap

	k=0;
	for(i=0;i<numspots;i++)
	{
		p[3*i]=RND*star->r*0.1;   //radius -prefer small spots to start 
//		p[3*i]=RND*star->r;   //radius -original version
		if(!FIXTHETAS)
			p[3*i+1]=acos(2.0*RND-1.0);	 //theta
		p[3*i+2]=RND*PIt2; //phi
		
		rspot=p[3*i];
		alpha[i]=asin(rspot/star->r);  
		theta=p[3*i+1];
		phi=p[3*i+2];
		while(theta>PI)
			theta-=PIt2;
		while(theta<(-PI))
			theta+=PIt2;
		if(theta<0)
		{
			theta=(-theta);
			phi+=PI;
		}
		while(phi>PIt2)
			phi-=PIt2;
		while(phi<0)
			phi+=PIt2;
		posunit[i][0]=sin(theta)*cos(phi); 
		posunit[i][1]=sin(theta)*sin(phi);
		posunit[i][2]=cos(theta);
		for(j=0;j<i;j++)
		{
			ang=acos(posunit[i][0]*posunit[j][0]+posunit[i][1]*posunit[j][1]+posunit[i][2]*posunit[j][2]);
			if(ang<alpha[i]+alpha[j])
			{				//redo point if it overlaps
				i--;
				j=i+10;
				k++;
				if(k>20&&i>=0)
				{			//too many redo's, go back further
					i--;
					k=0;
				}
			}
		}
		if(j<i+5)
			k=0;
	}
	if(CALCBRIGHTNESSFACTOR)
		p[3*numspots]=1.0;
	else
		p[3*numspots]=1.2-RND*0.4;  //star brightness factor
}
double modelmaxlight(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS])
{
//	FILE *debugmml;
	int realprintvis,realnumplanets,i,j,k,nt,localmaxind[2];
	double period,t,dt,light[360],localmaxl[2],t0,t1,l0,l1;
//	debugmml=fopen("debugmml.txt","w");
	star->brightnessfactor=1.0;
	realnumplanets=numplanets;
	numplanets=0;
	period=PIt2/star->omega;
	nt=360;
	dt=period/nt;
	realprintvis=PRINTVIS;
	PRINTVIS=0;
	for(i=0;i<360;i++)
	{
		t=dt*((double)i);
		light[i]=lightness(t,star,planet,spot);
//		fprintf(debugmml,"%i %lf %lf\n",i,t,light[i]);
	}
	localmaxl[0]=localmaxl[1]=0.0;
	localmaxind[0]=localmaxind[1]=(-1);
	if(light[0]>light[1]&&light[0]>=light[nt-1])
	{
		localmaxl[0]=light[0];
		localmaxind[0]=0;
	}
	if(light[nt-1]>=light[nt-2]&&light[nt-1]>light[0])
	{
		if(light[nt-1]>localmaxl[0])
		{
			localmaxl[1]=localmaxl[0];
			localmaxind[1]=localmaxind[0];
			localmaxl[0]=light[nt-1];
			localmaxind[0]=nt-1;
		}
		else
		{
			localmaxl[1]=light[nt-1];
			localmaxind[1]=nt-1;
		}
	}
	for(i=1;i<nt-1;i++)
		if(light[i]>=light[i-1]&&light[i]>light[i+1])
			if(light[i]>localmaxl[1])
				if(light[i]>localmaxl[0])
				{
					localmaxl[1]=localmaxl[0];
					localmaxind[1]=localmaxind[0];
					localmaxl[0]=light[i];
					localmaxind[0]=i;
				}
				else
				{
					localmaxl[1]=light[i];
					localmaxind[1]=i;
				}
//	fprintf(debugmml,"\nlocal max:\n%i %lf\n%i %lf\n\n",localmaxind[0],localmaxl[0],localmaxind[1],localmaxl[1]);
	for(i=0;i<2&&localmaxind[i]>=0;i++)
	{
		if(localmaxind[i]==0)
			j=nt-1;
		else
			j=localmaxind[i]-1;
		if(localmaxind[i]==nt-1)
			k=0;
		else
			k=localmaxind[i]+1;
		if(light[j]>light[k])
		{
			t0=((double)(localmaxind[i]-1))*dt;
			t1=((double)localmaxind[i])*dt;
			l0=light[j];
			l1=localmaxl[i];
		}
		else
		{
			t0=((double)localmaxind[i])*dt;
			t1=((double)(localmaxind[i]+1))*dt;
			l0=localmaxl[i];
			l1=light[k];
		}
		for(j=0;j<10;j++)
		{
//			fprintf(debugmml,"(%0.10lf %0.10lf) - (%0.10lf %0.10lf)\n",t0,l0,t1,l1);
			t=(t0+t1)/2.0;
			if(l0<l1)
			{
				t0=(t0+t1)/2.0;
				l0=lightness(t0,star,planet,spot);
			}
			else
			{
				t1=(t0+t1)/2.0;
				l1=lightness(t,star,planet,spot);
			}
		}
		if(l0>l1)
			localmaxl[i]=l0;
		else
			localmaxl[i]=l1;
//		fprintf(debugmml,"localmaxl[%i]= %lf\n\n",i,localmaxl[i]);
	}
	numplanets=realnumplanets;
	PRINTVIS=realprintvis;
//	fclose(debugmml);
	if(localmaxl[0]==0.0&&localmaxl[1]==0.0&&localmaxind[0]==(-1)&&localmaxind[1]==(-1))
		return light[0];	//all light[] values are the same
	if(localmaxl[0]>localmaxl[1])
		return localmaxl[0];
	else
	{
//		fprintf(outerr,"used localmax 1\n");
		return localmaxl[1];
	}
}
char setspots(double p[],spotdata spot[MAXSPOTS],stardata *star,planetdata planet[MAXPLANETS])
{                   //sets spots for call to lightness, also sets star brightness factor
	int i,j;
	double ang,posunit[MAXSPOTS][3];

	for(i=0;i<numspots;i++)
	{
		spot[i].fracdarkness=spotfracdark;
		spot[i].r=p[3*i];
		if(spot[i].r<0||spot[i].r>star->r)
			return 0;  //bad params
		spot[i].thetarel=p[3*i+1];
		if(!(spot[i].thetarel>=0&&spot[i].thetarel<=PI))
			return 0;	//bad params
		spot[i].phirel=p[3*i+2];
		if(!(spot[i].phirel>=0&&spot[i].phirel<=PIt2))
			return 0;	//bad params
		while(spot[i].thetarel>PI)
			spot[i].thetarel-=PIt2;
		while(spot[i].thetarel<(-PI))
			spot[i].thetarel+=PIt2;
		if(spot[i].thetarel<0)
		{
			spot[i].thetarel=(-spot[i].thetarel);
			spot[i].phirel+=PI;
		}
		while(spot[i].phirel>PIt2)
			spot[i].phirel-=PIt2;
		while(spot[i].phirel<0)
			spot[i].phirel+=PIt2;
		spot[i].alpha=asin(spot[i].r/star->r);
		spot[i].flatdcen=sqrt(star->r*star->r-spot[i].r*spot[i].r);
		spot[i].area=PI*spot[i].r*spot[i].r;
		spot[i].poscoeff[0]=(-sin(spot[i].thetarel)*sin(star->theta));
		spot[i].poscoeff[1]=cos(spot[i].thetarel)*cos(star->theta);
		spot[i].poscoeff[2]=sin(spot[i].thetarel)*cos(star->theta);
		spot[i].poscoeff[3]=cos(spot[i].thetarel)*sin(star->theta);
		spot[i].poscoeff[4]=sin(spot[i].thetarel);
		
		posunit[i][0]=sin(spot[i].thetarel)*cos(spot[i].phirel); 
		posunit[i][1]=sin(spot[i].thetarel)*sin(spot[i].phirel);
		posunit[i][2]=cos(spot[i].thetarel);
		for(j=0;j<i;j++)
		{
			ang=acos(posunit[i][0]*posunit[j][0]+posunit[i][1]*posunit[j][1]+posunit[i][2]*posunit[j][2]);
			if(ang<spot[i].alpha+spot[j].alpha)
				return 0;   //bad params
		}
	}
	if(CALCBRIGHTNESSFACTOR)
		p[3*numspots]=PI/modelmaxlight(star,planet,spot);
	star->brightnessfactor=p[3*numspots];
	return 1;
}
void combineangles(double theta0,double phi0,double theta1,double phi1,double z,double thphi2[2])
{
	double cth0,sth0,cth1,sth1,sth1cdf,th1pp,sth1pp,th2pp,sth2pp,cth2pp,cf2pp,sf2pp,sth2p,f2p,df;

	if(theta0==theta1&&phi0==phi1)
	{
		thphi2[0]=theta0;
		thphi2[1]=phi0;
	}
	else if(FIXTHETAS)
	{
		df=phi1-phi0;
		if(df<-PI)
			df+=PIt2;
		if(df>PI)
			df-=PIt2;
		thphi2[0]=theta0;
		thphi2[1]=phi1-z*df;
	}
	else
	{
		cth0=cos(theta0);
		sth0=sin(theta0);
		cth1=cos(theta1);
		sth1=sin(theta1);
		sth1cdf=sth1*cos(phi1-phi0);
		th1pp=acos(sth0*sth1cdf+cth0*cth1);
		sth1pp=sin(th1pp);
		th2pp=th1pp*z;
		sth2pp=sin(th2pp);
		cth2pp=cos(th2pp);
		cf2pp=(cth0*sth1cdf-sth0*cth1)/sth1pp;
		sf2pp=sth1*sin(phi1-phi0)/sth1pp;
		thphi2[0]=acos(cth0*cth2pp-sth0*sth2pp*cf2pp);
		sth2p=sin(thphi2[0]);
		f2p=acos((cth0*sth2pp*cf2pp+sth0*cth2pp)/sth2p);
		if(sth2pp*sf2pp/sth2p<0)
			f2p=PIt2-f2p;
		thphi2[1]=f2p+phi0;
	}

	if(thphi2[0]<0) thphi2[0]+=PIt2;
	else if(thphi2[0]>PIt2) thphi2[0]-=PIt2;
	if(thphi2[1]<0) thphi2[1]+=PIt2;
	else if(thphi2[1]>PIt2) thphi2[1]-=PIt2;
	
}
double combiners(double r0,double theta0,double r1,double theta1,double theta2,double z,stardata *star)
{
#	if COMBINERSINTHETA
		double rtheta0,rtheta1,rtheta2;
		double psize0,psize1,psize2,r2;
		rtheta0=theta0+star->theta;
		rtheta1=theta1+star->theta;
		rtheta2=theta2+star->theta;
		if(r0>star->r*sin(rtheta0)||r1>star->r*sin(theta1))
			r2=r0+z*(r1-r0);
		else
		{
			psize0=r0*r0*sin(rtheta0);
			psize1=r1*r1*sin(rtheta1);
			psize2=psize0+z*(psize1-psize0);
			if(psize2<0)
				r2=(-1);  //setspots doesn't notice -nan as not a good param, so give it -1
			else
				r2=sqrt(psize2/sin(rtheta2));
		}

		return r2;
#	else
		return r0+z*(r1-r0);
#	endif
}
#if COMBINEONESPOT
void combineonespot(double p0[],double p1[],double p2[],double z)
{
	int i,j;
	double thphi[2];

	if(numseeded==0||FIXSEEDEDONLYPHI)
		i=rand()%numspots;
	else
		i=numseeded+rand()%(numspots-numseeded);
	for(j=0;j<numspots;j++)
		if(j!=i)
		{
			p2[j*3]=p1[j*3];
			p2[j*3+1]=p1[j*3+1];
			p2[j*3+2]=p1[j*3+2];
		}
	p2[i*3]=p0[i*3]+z*(p1[i*3]-p0[i*3]);
	combineangles(p0[i*3+1],p0[i*3+2],p1[i*3+1],p1[i*3+2],z,thphi);
	p2[i*3+1]=thphi[0];
	p2[i*3+2]=thphi[1];
	if(p2[i*3+1]<0)
	{
		p2[i*3+1]=(-p2[i*3+1]);
		p2[i*3+2]+=PI;
	}
	if(p2[i*3+1]>PI)
	{
		p2[i*3+1]=PIt2-p2[i*3+1];
		p2[i*3+2]+=PI;
	}
	while(p2[i*3+2]<0)
		p2[i*3+2]+=PIt2;
	while(p2[i*3+2]>PIt2)
		p2[i*3+2]-=PIt2;
}
#endif
void combineparam(double p0[],double p1[],double p2[],double z,stardata *star)
{
	int i;
	double thphi[2];

#	if COMBINEONESPOT
		if(CALCBRIGHTNESSFACTOR||rand()%(numspots+1)>0)
			combineonespot(p0,p1,p2,z);
		else
		{
			for(i=0;i<numspots*3;i++)
				p2[i]=p1[i];
			p2[numspots*3]=p0[numspots*3]+z*(p1[numspots*3]-p0[numspots*3]);
		}
#	else
		for(i=0;i<numspots;i++)
		{
			combineangles(p0[i*3+1],p0[i*3+2],p1[i*3+1],p1[i*3+2],z,thphi);
			p2[i*3+1]=thphi[0];
			p2[i*3+2]=thphi[1];
			p2[i*3]=combiners(p0[i*3],p0[i*3+1],p1[i*3],p1[i*3+1],p2[i*3+1],z,star);
			if(p2[i*3+1]<0)
			{
				p2[i*3+1]=(-p2[i*3+1]);
				p2[i*3+2]+=PI;
			}
			if(p2[i*3+1]>PI)
			{
				p2[i*3+1]=PIt2-p2[i*3+1];
				p2[i*3+2]+=PI;
			}
			while(p2[i*3+2]<0)
				p2[i*3+2]+=PIt2;
			while(p2[i*3+2]>PIt2)
				p2[i*3+2]-=PIt2;
		}
		if(!CALCBRIGHTNESSFACTOR)
			p2[numspots*3]=p0[numspots*3]+z*(p1[numspots*3]-p0[numspots*3]);
#	endif
}
int orderspots(double param[])
{
	int j,k,cts;
	int of0,of1;
	double tmp;
	double thetap0,thetap1,dtheta,dphi;

	if(ORDERSPOTSMETHOD<2)
	{
		if(ORDERSPOTSMETHOD==1)
		{
			of0=2;	//order by longitude (phi)
			of1=1;
		}
		else
		{
			of0=1;	//order by latitude (theta)
			of1=2;
		}

		cts=0;
		j=numseeded;  
		while(j<numspots-1)
		{
			if(param[j*3+of0]<param[j*3+of0+3]-ORDERSPOTNEARNESS)
				j++;
			else if((param[j*3+of0]<param[j*3+of0+3]+ORDERSPOTNEARNESS)&&(param[j*3+of1]<param[j*3+of1+3]))
				j++;
			else
			{
				cts=1;
				for(k=0;k<3;k++)
				{
					tmp=param[j*3+k];
					param[j*3+k]=param[j*3+3+k];
					param[j*3+3+k]=tmp;
				}
				if(j>numseeded)
					j--;
				else
					j++;
			}
		}	
	}
	else
	{	//hybrid ordering method
		dphi=PI/18.0;
		cts=0;
		j=numseeded;
		while(j<numspots-1)
		{
			thetap0=param[j*3+1];
			if(thetap0>PIo2)	
				thetap0=PI-thetap0;
			thetap1=param[j*3+4];
			if(thetap1>PIo2)
				thetap1=PI-thetap1;
			dtheta=PI/18.0+0.833333*thetap0;	// pi/18 at theta=0 to pi/3 at theta=3
			if(thetap0<thetap1-dtheta)
				j++;
			else if(thetap0>thetap1+dtheta)
			{
				cts=1;
				for(k=0;k<3;k++)
				{
					tmp=param[j*3+k];
					param[j*3+k]=param[j*3+3+k];
					param[j*3+3+k]=tmp;
				}
				if(j>numseeded)
					j--;
				else
					j++;
			}
			else	
			{	//same theta band
				if(param[j*3+2]<param[j*3+5]-dphi)
					j++;
				else if(param[j*3+2]>param[j*3+5]+dphi)
				{
					cts=1;
					for(k=0;k<3;k++)
					{
						tmp=param[j*3+k];
						param[j*3+k]=param[j*3+3+k];
						param[j*3+3+k]=tmp;
					}
					if(j>numseeded)
						j--;
					else
						j++;
				}
				else
				{	//same theta and phi band
					if(thetap0<=thetap1)
						j++;
					else
					{
						cts=1;
						for(k=0;k<3;k++)
						{
							tmp=param[j*3+k];
							param[j*3+k]=param[j*3+3+k];
							param[j*3+3+k]=tmp;
						}
						if(j>numseeded)
							j--;
						else
							j++;
					}
				}
			}
		}	
	}

	return cts;	
}
char rndvarparam(double param[],double rndparam[],double sigmaspotradius,double sigmaspotpos,stardata *star)
{
	int i;
	double sth0,cth0,thr,fr,p1ppx,p1ppy,p1ppz,sthf;

	for(i=0;i<numspots;i++)
	{
		rndparam[i*3]=normaldist(param[i*3],sigmaspotradius); 
		if(rndparam[i*3]<0)
			return 0;
		if(FIXTHETAS)
		{
			rndparam[i*3+1]=param[i*3+1];
			rndparam[i*3+2]=normaldist(param[i*3+2],sigmaspotpos);
		}
		else if(FIXSEEDEDONLYPHI)
		{
			if(i<numseeded)
			{
				rndparam[i*3+1]=normaldist(param[i*3+1],sigmaspotpos);
				if(FIXSEEDEDONLYPHI>1)  
					rndparam[i*3+2]=normaldist(param[i*3+2],sigmaspotpos);  //fix nothing
				else
					rndparam[i*3+2]=param[i*3+2];  //fix only phi
			}
			else
			{
				rndparam[i*3]=param[i*3];
				rndparam[i*3+1]=param[i*3+1];
				rndparam[i*3+2]=param[i*3+2];
			}
		}
		else
		{
			sth0=sin(param[i*3+1]);
			cth0=cos(param[i*3+1]);
			thr=normaldist(0.0,sigmaspotpos);
			fr=RND*PIt2;
			p1ppx=sin(thr)*cos(fr);
			p1ppy=sin(thr)*sin(fr);
			p1ppz=cos(thr);
			rndparam[i*3+1]=acos(cth0*p1ppz-sth0*p1ppx);
			sthf=sin(rndparam[i*3+1]);
			rndparam[i*3+2]=acos((cth0*p1ppx+sth0*p1ppz)/sthf);
			if(p1ppy/sthf<0)
				rndparam[i*3+2]=PIt2-rndparam[i*3+2];
			rndparam[i*3+2]+=param[i*3+2];
		}
		while(rndparam[i*3+2]<0.0)
			rndparam[i*3+2]+=PIt2;
		while(rndparam[i*3+2]>PIt2)
			rndparam[i*3+2]-=PIt2;
	}
	rndparam[numspots*3]=normaldist(param[numspots*3],0.002);  //need to input sigma brightness correction -doesn't matter when CALCBRIGHTNESSFACTOR
	return 1;
}
void metrohast(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],double readparam[],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,int maxsteps,char rootname[64])
{
	double torig;
	char goodparams,filename[128];
	int i,j,msi,nparam;
	double *param,chisq,*potparam,potchisq,dchisq,r;
	FILE *bestparam,*bestlc;

	sprintf(filename,"%s_mhparambest.txt",rootname);
	bestparam=fopen(filename,"w");
	sprintf(filename,"%s_mhlcbest.txt",rootname);
	bestlc=fopen(filename,"w");

	nparam=3*numspots+1;
	param=(double *)malloc(nparam*sizeof(double));
	potparam=(double *)malloc(nparam*sizeof(double));
	if(param==NULL||potparam==NULL)
		printf("malloc error\n");
	for(i=0;i<nparam;i++)
		param[i]=readparam[i];
	goodparams=setspots(param,spot,star,planet);
	if(!goodparams)
		printf("mh error\n");
	chisq=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);

	printf("initial chisq: %lf\n",chisq);
	for(msi=0;msi<maxsteps;msi++)
	{
		printf(".");
		fflush(stdout);
		goodparams=rndvarparam(param,potparam,sigmaradius,sigmaangle,star);
		if(goodparams)
		{
			goodparams=setspots(potparam,spot,star,planet);
			if(goodparams)
			{
				potchisq=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				if(potchisq>=0)
				{
					printf("\n%lf\n",potchisq);
					dchisq=potchisq-chisq;
					j=0;
					if(dchisq<0)
						j=1;
					else
					{
						r=exp(-dchisq/2.0);
						if(RND<r)
							j=1;
					}
					if(j)	//accept
					{
						for(i=0;i<nparam;i++)
							param[i]=potparam[i];
						chisq=potchisq;
						printf("\n%lf   (%i)\n",chisq,msi);
					}
				}
			}
		}
	}

	for(j=0;j<nparam;j++)
		fprintf(bestparam,"%lf\n",param[j]);
	fprintf(bestparam,"%lf    chi squared\n",chisq);
	(void)setspots(param,spot,star,planet);
	for(j=0;j<numspots;j++)
		fprintf(bestparam,"%lf    radius\n%lf    theta\n%lf    phi\n",spot[j].r,spot[j].thetarel*180.0/PI,spot[j].phirel*180/PI);

#	if ANYPRINTVIS&&ALWAYSPRINTVIS
		PRINTVIS=WHICHPRINTVIS;
		sprintf(filename,"%s_vis.txt",rootname);
		outv=fopen(filename,"w");
#	endif

	for(i=0;i<lcn;i++)
	{
		torig=lctime[i]/86400.0+planet[0].lct0;
		r=lightness(lctime[i],star,planet,spot);
#		if UNNORMALIZEOUTPUT
			fprintf(bestlc,"%0.9lf %0.6lf %0.6lf %0.6lf",torig,lclight[i]/lclightnorm,lcuncertainty[i]/lclightnorm,r/lclightnorm);
#			if PRINTEACHSPOT
				for(j=0;j<numspots;j++)
					fprintf(bestlc," %0.6lf",spotreport[j]/lclightnorm);
#			endif
			fprintf(bestlc,"\n");
#		else
			fprintf(bestlc,"%0.9lf %0.6lf %0.6lf %0.6lf",torig,lclight[i],lcuncertainty[i],r);
#			if PRINTEACHSPOT
				for(j=0;j<numspots;j++)
					fprintf(bestlc," %0.6lf",spotreport[j]);
#			endif
			fprintf(bestlc,"\n");
#		endif
#		if ANYPRINTVIS
			if(PRINTVIS&&PRINTVIS!=1)
				fprintf(outv,"z\n%lf %lf\n",lclight[i],r);
#		endif
	}

	fclose(bestparam);
	fclose(bestlc);
	free((void *)param);
	free((void *)potparam);

}
double plotdataopt(int varyspot,int wntv,double hrange,double *param,stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm)
{
	int TIMESTOBIF=10;  //how many times to split the segment
	char goodparams;
	double best,bestp[2],chisqt[11],startp;
	int i,j,vind,besti=0;

	best=TOOBIG-1;
	vind=3*varyspot+wntv;
	startp=param[vind];
	for(i=(-5);i<=5;i++)
	{
		param[vind]=startp+((double)i/5.0)*hrange;
		goodparams=setspots(param,spot,star,planet);
		if(goodparams)
			chisqt[i+5]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
		else
			chisqt[i+5]=TOOBIG;
		if(chisqt[i+5]<best)
		{
			besti=i;
			best=chisqt[i+5];
		}
	}

	if(besti==(-5))
	{
		bestp[0]=startp-hrange;
		bestp[1]=startp-0.8*hrange;
	}
	else if(besti==5)
	{
	
		bestp[0]=startp+0.8*hrange;
		bestp[1]=startp+hrange;
		chisqt[0]=chisqt[9];
		chisqt[1]=chisqt[10];
	}
	else if(chisqt[besti+4]<chisqt[besti+6])
	{
		bestp[0]=startp+((double)(besti-1)/5.0)*hrange;
		bestp[1]=startp+((double)besti/5.0)*hrange;
		chisqt[0]=chisqt[besti+4];
		chisqt[1]=chisqt[besti+5];
	}
	else
	{
		bestp[0]=startp+((double)besti/5.0)*hrange;
		bestp[1]=startp+((double)(besti+1)/5.0)*hrange;
		chisqt[0]=chisqt[besti+5];
		chisqt[1]=chisqt[besti+6];
	}

	for(i=0;i<TIMESTOBIF;i++)
	{
		if(chisqt[0]<chisqt[1])
			j=1;
		else
			j=0;
		bestp[j]=0.5*(bestp[0]+bestp[1]);
		param[vind]=bestp[j];
		goodparams=setspots(param,spot,star,planet);
		if(goodparams)
			chisqt[j]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
		else
			chisqt[j]=TOOBIG;
	}

	if(chisqt[0]<chisqt[1])
	{
		param[vind]=bestp[0];
		return chisqt[0];
	}
	else
	{
		param[vind]=bestp[1];
		return chisqt[1];
	}
}
void plotdatacolor(char colorscheme, double x, int *r, int *g,int *b)
{
	int red=0,green=0,blue=0;
	double col[4];

	if(colorscheme==0)
	{
		col[0]=100.0;	
		col[1]=600.0;
		col[2]=1200.0;
		col[3]=5000.0;

		red=255-((int)(x*255.0/col[0]));
		if(red<0)
			red=0;
		if(x<col[1])
			green=255-(int)(255.0*(col[1]-x)/(col[1]-col[0]));
		else
			green=255-(int)(255.0*(x-col[1])/(col[2]-col[1]));
		if(green<0)
			green=0;
		blue=255-(int)(255.0*(col[3]-x)/(col[3]-col[2]));
		if(blue<0)
			blue=0;
		if(blue>255)
			blue=255;
	}
	else if(colorscheme==1)
	{
		col[0]=2;
		col[1]=30;

		if(x*col[0]<255)
			red=(int)(255.0-x*col[0]);
		else
			green=(int)(col[1]*log10(col[0]*x-250));
		if(red<0)
			red=0;
		if(red>255)
			red=255;
		if(green<0)
			green=0;
		if(green>255)
			green=255;
	}

	*r=red;
	*g=green;
	*b=blue;
}
void plotdata(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,int varyspot,double readparam[],char rootname[64])
{
	char goodparams,filename[128],innotout;
	char wtv[2],wntv,optimize;  //which parameter to vary,which not to vary (the other one), whether to optimize the unvaried parameter
	int i,j,n,nparam,red,green,blue;
	int hgrid[3],w,h;  //number of points to vary to, half grid
	double hrange[3];  //half range of variation
	double *param,chisq; //param[0-1][chain number][parameter number]
	FILE *out,*outppm,*in;

	innotout=1;  //1 for reading previous run's data and changing the colors

	sprintf(filename,"%s_%i_plotdata.txt",rootname,varyspot);
	if(innotout)
		in=fopen(filename,"r");
	else
		out=fopen(filename,"w");
	sprintf(filename,"%s_%i.ppm",rootname,varyspot);
	outppm=fopen(filename,"w");

	nparam=3*numspots+1; //r, theta, phi of spots, then unoccluded star brightness factor
	param=(double *)malloc(nparam*sizeof(double));
		if(param==NULL)
			printf("malloc error\n");

	for(i=0;i<nparam;i++)
		param[i]=readparam[i];
	goodparams=setspots(param,spot,star,planet);
	if(!goodparams)
	{
		printf("bad seed\n");
		exit(0);
	}
	chisq=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
	if(chisq<0)
	{
		printf("bad seed\n");
		exit(0);
	}

	wtv[0]=1;  //what parameters to vary and what not to, 0=raduis, 1=theta, 2=phi
	wtv[1]=0;
	wntv=2;
	optimize=0;  //1->optimize unvaried parameter, 0->use seed value

	hgrid[0]=250;		//radius			need range for optimized parameters as well
	hrange[0]=0.1;
	hgrid[1]=250;	//theta
	hrange[1]=0.5;
	hgrid[2]=0;	//phi
	hrange[2]=0.0;

	h=2*hgrid[wtv[0]]+1;
	w=2*hgrid[wtv[1]]+1;
	fprintf(outppm,"P3\n%i %i\n255\n",w,h);

	n=0;
	goodparams=1;
	for(i=-hgrid[wtv[0]];i<=hgrid[wtv[0]];i++)
		for(j=-hgrid[wtv[1]];j<=hgrid[wtv[1]];j++)
		{ 
			param[varyspot*3+wtv[0]]=readparam[varyspot*3+wtv[0]]+((double)i)*((double)hrange[wtv[0]])/hgrid[wtv[0]];
			param[varyspot*3+wtv[1]]=readparam[varyspot*3+wtv[1]]+((double)j)*((double)hrange[wtv[1]])/hgrid[wtv[1]];
			param[varyspot*3+wntv]=readparam[varyspot*3+wntv];

			if(innotout)
				fscanf(in,"%lf %lf %lf %lf\n",param+varyspot*3,param+varyspot*3+1,param+varyspot*3+2,&chisq);
			else
			{
				if(optimize)
					chisq=plotdataopt(varyspot,wntv,hrange[wntv],param,star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm);
				else
				{
					goodparams=setspots(param,spot,star,planet);
					if(goodparams)
						chisq=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				}
				fprintf(out,"%lf %lf %lf %lf\n",param[varyspot*3],param[varyspot*3+1],param[varyspot*3+2],chisq);
			}
			if(goodparams)
			{
				plotdatacolor(1,chisq,&red,&green,&blue);
				fprintf(outppm,"%i %i %i\n",red,green,blue);
				printf(".");
			}
			else
			{
				printf("x");
				red=255;
				blue=255;
				green=255;
				fprintf(outppm,"%i %i %i\n",red,green,blue);
			}
			n++;
			if(n==25)
			{
				printf(" %i %i\n",i,j);
				n=0;
			}
		}
		
	if(innotout)
		fclose(in);
	else
		fclose(out);
	free((void *)param);
}
void mcmc(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,double ascale,int npop,long int maxstepsortime,char rootname[64],char seeded,double readparam[],char seedfilename[64])
{
	char goodparams,filename[128];
	int i,j,k,r,msi,nparam,curstep,potstep,*updated,*naccepted,maxsteps;
	long int memused,maxtime,time0,time1,avgtime;
	double **param[2],*chisq[2]; //param[0-1][chain number][parameter number]
	double sqrtascale,oosqrtascale,smooascale,z,alpha,rn;
	double bestchisq,*bestparam,torig;
	unsigned int whichspot;
	FILE *parambest,*bestlc,*tracker,*finalparam,*in;
#	if MCMCTRACKMEM
		int memtrackind,*memtracki;
		double *memtrackd;
#	endif
#	if DEBUGMCMC
		int d;
#	endif

	nparam=3*numspots+1; //r, theta, phi of spots, then unoccluded star brightness factor

#	if MCMCTRACKMEM
		memtracki=(int *)malloc(3*MCMCTRACKMEM*sizeof(int));
		memtrackd=(double *)malloc((nparam+1)*MCMCTRACKMEM*sizeof(double));
		if(memtracki==NULL||memtrackd==NULL)
		{
			printf("malloc error for tracking buffer\n");
			exit(0);
		}
		memtrackind=0;
#	endif
#	if !QUIET && MCMCTRACKMEM
		printf("memory used for tracker buffer: %lu\n",3*MCMCTRACKMEM*sizeof(int)+(nparam+1)*MCMCTRACKMEM*sizeof(double));
#	endif

	sprintf(filename,"%s_parambest.txt",rootname);
	parambest=fopen(filename,"w");
	if(parambest==NULL)
		printf("error opening parambest file: %s\n",filename);
	sprintf(filename,"%s_lcbest.txt",rootname);
	bestlc=fopen(filename,"w");
	if(bestlc==NULL)
		printf("error opening bestlc file: %s\n",filename);
	sprintf(filename,"%s_mcmc.txt",rootname);
	tracker=fopen(filename,"w");
	if(tracker==NULL)
		printf("error opening tracker file: %s\n",filename);
	sprintf(filename,"%s_finalparam.txt",rootname);
	finalparam=fopen(filename,"w");
	if(finalparam==NULL)
		printf("error opening finalparam file: %s\n",filename);

	bestchisq=10000000000;

	if(seeded==2)
	{
		in=fopen(seedfilename,"r");
		if(in==NULL)
		{
			printf("no seed param file\n");
			exit(0);
		}
		fscanf(in,"%i\n%i\n",&i,&j);
		if(i!=numspots)
		{
#			if !QUIET
				printf("using numspots= %i from file\n",i);
#			endif
			numspots=i;
		}
		if(j!=npop)
		{
#			if !QUIET
				printf("using npop= %i from file\n",j);
#			endif
			npop=j;
		}
	}

	sqrtascale=sqrt(ascale);
	oosqrtascale=1.0/sqrtascale;
	smooascale=sqrtascale-oosqrtascale;
	updated=(int *)malloc(npop*sizeof(int));
	naccepted=(int *)malloc(npop*sizeof(int));
	chisq[0]=(double *)malloc(npop*sizeof(double));
	chisq[1]=(double *)malloc(npop*sizeof(double));
	param[0]=(double **)malloc(npop*sizeof(double *));
	param[1]=(double **)malloc(npop*sizeof(double *));
	memused=2*npop*sizeof(int)+2*npop*sizeof(double)+2*npop*sizeof(double *);
	for(i=0;i<npop;i++)
	{
		naccepted[i]=0;
		updated[i]=0;
		for(j=0;j<2;j++)
		{
			param[j][i]=(double *)malloc(nparam*sizeof(double));
			if(param[j][i]==NULL)
				printf("malloc error\n");
			memused+=nparam*sizeof(double);
		}
	}
#	if !QUIET
		printf("memory used for mcmc: %li\n",memused);
#	endif
	bestparam=(double *)malloc(nparam*sizeof(double));

#	if !QUIET
		if(FIXTHETAS)
			printf("using fixed values for spot thetas\n");
#	endif
	if(!seeded)	//not seeded
	{
		if(FIXTHETAS)
			for(i=0;i<npop;i++)
				for(j=0;j<numspots;j++)
					param[0][i][j*3+1]=readparam[j];
		for(i=0;i<npop;i++)
		{
			setrandomparam(param[0][i],star);
			(void)setspots(param[0][i],spot,star,planet);
			chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
			while(chisq[0][i]==(-1))
			{
				setrandomparam(param[0][i],star);
				(void)setspots(param[0][i],spot,star,planet);
				chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
			}
		}
	}
	else if(seeded==1)	//seeded
	{
		for(i=0;i<nparam;i++)
			param[0][0][i]=readparam[i];
		goodparams=setspots(param[0][0],spot,star,planet);
		if(!goodparams)
		{
			printf("bad seed\n");
			exit(0);
		}
		chisq[0][0]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
		if(chisq[0][0]<0)
		{
			printf("bad seed\n");
			exit(0);
		}
		for(i=1;i<npop;i++)
		{
			goodparams=0;
			while(!goodparams)
			{
				goodparams=rndvarparam(param[0][0],param[0][i],sigmaradius,sigmaangle,star);
				if(goodparams)
				{
					goodparams=setspots(param[0][i],spot,star,planet);
					if(goodparams)
					{
						chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
						if(chisq[0][i]>=0)
							goodparams=1;
						else
							goodparams=0;
					}
				}
			}
		}
	}
	else if(seeded==2)	//totally seeded
	{
		for(j=0;j<npop;j++)
		{
			for(i=0;i<nparam;i++)
				fscanf(in,"%lf\n",&param[0][j][i]);
			goodparams=setspots(param[0][j],spot,star,planet);
			if(!goodparams)
			{
				printf("bad seed\n");
				exit(0);
			}
			chisq[0][j]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
			if(chisq[0][j]<0)
			{
				printf("bad seed\n");
				exit(0);
			}
		}

		fclose(in);
	}
	else	//partially seeded
	{
		for(i=0;i<npop;i++)
		{
			if(FIXSEEDEDONLYPHI)
			{
				setrandomparam(param[1][i],star);
				for(j=0;j<numseeded*3;j++)
					param[1][i][j]=readparam[j];
				(void)rndvarparam(param[1][i],param[0][i],sigmaradius,sigmaangle,star);
				goodparams=setspots(param[0][i],spot,star,planet);
				if(!goodparams)
					chisq[0][i]=(-1);
				else
					chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				while(chisq[0][i]==(-1))
				{
					setrandomparam(param[1][i],star);
					for(j=0;j<numseeded*3;j++)
						param[1][i][j]=readparam[j];
					(void)rndvarparam(param[1][i],param[0][i],sigmaradius,sigmaangle,star);
					goodparams=setspots(param[0][i],spot,star,planet);
					if(!goodparams)
						chisq[0][i]=(-1);
					else
						chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				}
			}
			else
			{
				setrandomparam(param[0][i],star);
				for(j=0;j<numseeded*3;j++)
					param[0][i][j]=readparam[j];
				goodparams=setspots(param[0][i],spot,star,planet);
				if(!goodparams)
					chisq[0][i]=(-1);
				else
					chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				while(chisq[0][i]==(-1))
				{
					setrandomparam(param[0][i],star);
					for(j=0;j<numseeded*3;j++)
						param[0][i][j]=readparam[j];
					goodparams=setspots(param[0][i],spot,star,planet);
					if(!goodparams)
						chisq[0][i]=(-1);
					else
						chisq[0][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				}
			}
		}
	}

	for(i=0;i<npop;i++)		//print out initial chains to _mcmc file
	{
		fprintf(tracker,"%i %i %i %lf",i,0,0,chisq[0][i]);
		for(j=0;j<nparam;j++)
			fprintf(tracker," %lf",param[0][i][j]);
		fprintf(tracker,"\n");
	}



	bestchisq=chisq[0][0];
	for(i=0;i<nparam;i++)
		bestparam[i]=param[0][0][i];
	for(j=1;j<npop;j++)
		if(chisq[0][j]<bestchisq)
		{
			bestchisq=chisq[0][j];
			for(i=0;i<nparam;i++)
				bestparam[i]=param[0][0][i];
		}

	if(maxstepsortime>0)
	{
		maxsteps=maxstepsortime;
		maxtime=0;
	}
	else
	{
		maxsteps=TOOBIG;
		maxtime=starttime-maxstepsortime;  //time limit is passed as negative
#		if !QUIETMCMC && !QUIET
			printf("start time: %li maxtime: %li   (%li)\n",starttime,maxtime,maxstepsortime*(-1));
#		endif
	}
	avgtime=0;

	for(msi=0;msi<maxsteps;msi++)	//the mcmc loop
	{
		if(maxtime||ALWAYSAVERAGETIME)
			time0=timecheck();
		j=0;
		for(i=0;i<npop;i++)
			j+=naccepted[i];
#		if !QUIETMCMC && !QUIET
			printf("%i  (%i)",msi,j);
			if(maxtime||ALWAYSAVERAGETIME)
				printf("  avg step time: %li\n",avgtime);
			else
				printf("\n");
#		endif
		if((msi<=100&&msi>=25&&msi%5==0)||(msi%1005==0&&msi>0))
		{
#			if !QUIETMCMC && !QUIET
				printf("reordering spots   ");
#			endif
			i=0;
			for(j=0;j<npop;j++)
				i+=orderspots(param[naccepted[j]%2][j]);
#			if !QUIETMCMC && !QUIET
				printf("%i\n",i);
#			endif
		}
		for(i=0;i<npop;i++)
		{
			curstep=naccepted[i]%2;
			potstep=1-curstep;
			r=(int)floor(RND*npop);
			while(r==i||r<0||r>=npop)
				r=(int)floor(RND*npop);
			z=RND*smooascale+oosqrtascale;
			z=z*z;
			if(r<i)
				k=(naccepted[r]+updated[r])%2;
			else
				k=naccepted[r]%2;
			combineparam(param[k][r],param[curstep][i],param[potstep][i],z,star);					
			goodparams=setspots(param[potstep][i],spot,star,planet);
			if(goodparams)
				chisq[potstep][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
			else
				chisq[potstep][i]=(-2);
			while(chisq[potstep][i]==(-1))
			{		//happens when it is findchisq's fault there was an error
				r=(int)floor(RND*npop);
				while(r==i||r<0||r>=npop)
					r=(int)floor(RND*npop);
				z=RND*smooascale+oosqrtascale;
				z=z*z;
				if(r<i)
					k=(naccepted[r]+updated[r])%2;
				else
					k=naccepted[r]%2;
				combineparam(param[k][r],param[curstep][i],param[potstep][i],z,star);					
				goodparams=setspots(param[potstep][i],spot,star,planet);
				if(goodparams)
					chisq[potstep][i]=findchisq(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
				else
					chisq[potstep][i]=(-2);
			}
#			if DEBUGMCMC
				if(i==5)
				{
					fprintf(dbmcmc,"%i %lf ",msi,z);
					for(d=0;d<nparam;d++)
						fprintf(dbmcmc,"%lf ",param[curstep][i][d]);
					fprintf(dbmcmc,"%lf ",chisq[curstep][i]);
					for(d=0;d<nparam;d++)
						fprintf(dbmcmc,"%lf ",param[k][r][d]);
					fprintf(dbmcmc,"%lf ",chisq[k][r]);
					for(d=0;d<nparam;d++)
						fprintf(dbmcmc,"%lf ",param[potstep][i][d]);
					fprintf(dbmcmc,"%lf ",chisq[potstep][i]);
				}
#			endif
			if(chisq[potstep][i]>=0)
			{
				alpha=pow(z,nparam-1)*exp(-0.5*(chisq[potstep][i]-chisq[curstep][i]));
				rn=RND;
			}
			else
			{
				alpha=0;	//make it reject
				rn=1;
			}

#			if DEBUGMCMC
				if(i==5)
				{
					if(alpha>=rn)
						fprintf(dbmcmc,"1\n");
					else
						fprintf(dbmcmc,"0\n");
				}
#			endif

			if(alpha>=rn)
			{	//accept step
				naccepted[i]++;
				updated[i]=1;
#				if !MCMCTRACKMEM
					fprintf(tracker,"%i %i %i %lf",i,naccepted[i],msi+1,chisq[potstep][i]);
					for(j=0;j<nparam;j++)
						fprintf(tracker," %lf",param[potstep][i][j]);
					fprintf(tracker,"\n");
#				else
					memtracki[memtrackind*3]=i;
					memtracki[memtrackind*3+1]=naccepted[i];
					memtracki[memtrackind*3+2]=msi+1;
					memtrackd[memtrackind*(nparam+1)]=chisq[potstep][i];
					for(j=0;j<nparam;j++)
						memtrackd[memtrackind*(nparam+1)+1+j]=param[potstep][i][j];
					memtrackind++;
					if(memtrackind==MCMCTRACKMEM)
					{
						for(memtrackind=0;memtrackind<MCMCTRACKMEM;memtrackind++)
						{
							fprintf(tracker,"%i %i %i",memtracki[memtrackind*3],memtracki[memtrackind*3+1],memtracki[memtrackind*3+2]);
							for(j=0;j<nparam+1;j++)
								fprintf(tracker," %lf",memtrackd[memtrackind*(nparam+1)+j]);
							fprintf(tracker,"\n");
						}
						memtrackind=0;
					}
#				endif
				if(chisq[potstep][i]<bestchisq)
				{
					bestchisq=chisq[potstep][i];
					for(j=0;j<nparam;j++)
						bestparam[j]=param[potstep][i][j];
#					if !QUIETMCMC && !QUIET
						printf("%i %i %lf\n",i,naccepted[i],bestchisq);
#					endif
				}	
			}
		}
		for(i=0;i<npop;i++)
			updated[i]=0;
		if(maxtime||ALWAYSAVERAGETIME)
		{
			time1=timecheck();
			avgtime=(avgtime*msi+time1-time0)/(msi+1);
			if(maxtime&&time1+2*avgtime+300>maxtime)
				maxsteps=0;			//make it stop, time is almost up
		}
	}			//end of main mcmc loop

#	if MCMCTRACKMEM
		k=memtrackind;
		for(memtrackind=0;memtrackind<k;memtrackind++)
		{
			fprintf(tracker,"%i %i %i",memtracki[memtrackind*3],memtracki[memtrackind*3+1],memtracki[memtrackind*3+2]);
			for(j=0;j<nparam+1;j++)
				fprintf(tracker," %lf",memtrackd[memtrackind*(nparam+1)+j]);
			fprintf(tracker,"\n");
		}
#	endif

	fprintf(finalparam,"%i\n%i\n",numspots,npop);
	k=0;
	for(i=0;i<npop;i++)
	{
		curstep=naccepted[i]%2;
		for(j=0;j<nparam;j++)
			fprintf(finalparam,"%0.12lf\n",param[curstep][i][j]);
	}

	for(j=0;j<nparam;j++)
		fprintf(parambest,"%0.12lf\n",bestparam[j]);
	fprintf(parambest,"%lf    chi squared\n",bestchisq);
	(void)setspots(bestparam,spot,star,planet);
	for(j=0;j<numspots;j++)
		fprintf(parambest,"%lf    radius\n%lf    theta\n%lf    phi\n",spot[j].r,spot[j].thetarel*180.0/PI,spot[j].phirel*180/PI);
	fprintf(parambest,"%lf     unoccluded brightness\n",star->brightnessfactor*star->area/lclightnorm);

#	if ANYPRINTVIS&&ALWAYSPRINTVIS
		PRINTVIS=WHICHPRINTVIS;
		sprintf(filename,"%s_vis.txt",rootname);
		outv=fopen(filename,"w");
#	endif

	for(i=0;i<lcn;i++)
	{
		torig=lctime[i]/86400.0+planet[0].lct0;
#		if PRINTWHICHSPOT
			for(j=0;j<16;j++)
				whichspots[j]=0;
#		endif
		z=lightness(lctime[i],star,planet,spot);
#		if UNNORMALIZEOUTPUT
			fprintf(bestlc,"%0.9lf %0.6lf %0.6lf %0.6lf",torig,lclight[i]/lclightnorm,lcuncertainty[i]/lclightnorm,z/lclightnorm);
#			if PRINTEACHSPOT
				for(j=0;j<numspots;j++)
					fprintf(bestlc," %0.6lf",spotreport[j]/lclightnorm);
#			endif
#			if PRINTWHICHSPOT
				whichspot=0;
				for(j=0;j<16;j++)
					whichspot+=whichspots[j]*ttt[j];
				fprintf(bestlc," %i",whichspot);
#			endif
			fprintf(bestlc,"\n");
#		else
			fprintf(bestlc,"%0.9lf %0.6lf %0.6lf %0.6lf",torig,lclight[i],lcuncertainty[i],z);
#			if(PRINTEACHSPOT)
				for(j=0;j<numspots;j++)
					fprintf(bestlc," %0.6lf",spotreport[j]);
#			endif
#			if PRINTWHICHSPOT
				whichspot=0;
				for(j=0;j<16;j++)
					whichspot+=whichspots[j]*ttt[j];
				fprintf(bestlc," %i",whichspot);
#			endif
			fprintf(bestlc,"\n");
#		endif
#		if ANYPRINTVIS
			if(PRINTVIS&&PRINTVIS!=1)
				fprintf(outv,"z\n%lf %lf\n",lclight[i],z);
#		endif
	}
	if(maxtime||ALWAYSAVERAGETIME)
		fprintf(outerr,"average step time: %li\n",avgtime);
	fclose(parambest);
	fclose(bestlc);
	fclose(tracker);
	fclose(finalparam);
	free((void *)updated);
	free((void *)naccepted);
	free((void *)chisq[0]);
	free((void *)chisq[1]);
	for(i=0;i<npop;i++)
		for(j=0;j<2;j++)
			free((void *)param[j][i]);
	free((void *)param[0]);
	free((void *)param[1]);

}
void lcgen(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,char lcoutfilename[128])
{
	int i,j,realnum;
	double light,torig,l1,l2;
	unsigned int whichspot;
	FILE *lcout;

	lcout=fopen(lcoutfilename,"w");
	for(i=0;i<lcn;i++)
	{
		debuglcni=i;
		errorflag=0;
		torig=lctime[i]/86400.0+planet[0].lct0;
#		if PRINTWHICHSPOT
			for(j=0;j<16;j++)
				whichspots[j]=0;
#		endif

		light=lightness(lctime[i],star,planet,spot);
#		if !UNNORMALIZEOUTPUT
			lclightnorm=1.0;
#		endif
		fprintf(lcout,"%0.9lf %0.6lf %0.6lf %0.6lf",torig,lclight[i]/lclightnorm,lcuncertainty[i]/lclightnorm,light/lclightnorm);
#		if PRINTEACHSPOT
			for(j=0;j<numspots;j++)
				fprintf(lcout," %0.6lf",spotreport[j]/lclightnorm);
#		endif
#		if PRINTWHICHSPOT
			whichspot=0;
			for(j=0;j<16;j++)
				whichspot+=whichspots[j]*ttt[j];
			fprintf(lcout," %i",whichspot);
#		endif
#		if PRINTPLANETSPOTOVERLAP
			realnum=numplanets;
			numplanets=0;
			l1=lightness(lctime[i],star,planet,spot);
			numplanets=realnum;
			realnum=numspots;
			numspots=0;
			l2=lightness(lctime[i],star,planet,spot);
			numspots=realnum;
			fprintf(lcout," %0.6lf",(light-l1-l2+star->area)/lclightnorm);
#		endif
		fprintf(lcout,"\n");
//			fprintf(lcout,"%0.9lf %0.6lf %0.6lf\n",torig,light/lclightnorm,lcuncertainty[i]/lclightnorm);	//for making test lightcurves
		if(errorflag)
		{
			fprintf(outerr,"error %i at datum %i (location 1)\n",errorflag,i);
			printf("lcgen: error\n");
			exit(0);
		}
#		if ANYPRINTVIS
			if(PRINTVIS)
			{
				if(PRINTVIS!=1)
					fprintf(outv,"z\n%lf %lf\n",lclight[i],light);
			}
#		endif
	}
	fclose(lcout);
}
void initialguess(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,char infilename[128])
{
	char outfilename[128];
	FILE *in,*out;
	int i,j;
	double tday,t,theta,phi;
	double r[3],xhat[3],zhat[3],dp;

	in=fopen(infilename,"r");

	j=4;
	for(i=0;i<128;i++)
		if(infilename[i]=='\0')
		{
			j=i;
			i=128;
		}
	if(infilename[j-4]=='.'&&infilename[j-3]=='t'&&infilename[j-2]=='x'&&infilename[j-1]=='t')
		infilename[j-4]='\0';
	sprintf(outfilename,"%s_guess.txt",infilename);
	out=fopen(outfilename,"w");

	zhat[0]=sin(star->theta);
	zhat[1]=0.0;
	zhat[2]=cos(star->theta);

	while(!feof(in))
	{
		fscanf(in,"%lf",&tday);
		t=(tday-planet[0].lct0)*86400.0;

		star->phi=star->omega*t;
		xhat[0]=cos(star->theta)*cos(star->phi);
		xhat[1]=sin(star->phi);
		xhat[2]=(-1.0)*sin(star->theta)*cos(star->phi);

		fprintf(out,"%0.9lf ",tday);
		for(i=0;i<numplanets;i++)
		{
			setplanetpos(t,planet,i);
			if(planet[i].x>0)	//planet is in front of star
			{
				if(sqrt(planet[i].y*planet[i].y+planet[i].z*planet[i].z)-planet[i].r<star->r)
				{
					r[0]=sqrt(star->rsq-planet[i].y*planet[i].y+planet[i].z*planet[i].z)/star->r;
					r[1]=planet[i].y/star->r;
					r[2]=planet[i].z/star->r;
					dp=r[0]*zhat[0]+r[1]*zhat[1]+r[2]*zhat[2];
					theta=acos(dp);
					dp=r[0]*xhat[0]+r[1]*xhat[1]+r[2]*xhat[2];
					phi=acos(dp/sin(theta));
					fprintf(out,"%0.9lf %0.9lf\n",theta,phi);
				}
				else
					fprintf(out,"planet center not over star\n");
			}
			else
				fprintf(out,"planet behind star\n");
		}
	}
	fclose(in);
	fclose(out);
}
int areacoveredupdate(int sn,char eorc,int np,int nq,vwpoint p[4],vwpoint q[2],double egamma[MAXSPOTS][6],double ccgamma[MAXSPOTS][4],double cegamma[MAXSPOTS][4],int ebeencovered[MAXSPOTS],int cbeencovered[MAXSPOTS],int enreg[MAXSPOTS],int ccregactive[MAXSPOTS][3],int ceregactive[MAXSPOTS][3])
{
	int i,j,k,nnr,cwr,ccwr;
	double cw[3],ccw[3];

	if(eorc=='e')
	{
		if(np==2)	//make into regions that don't cross over zero
		{
			if(p[0].gamma>p[1].gamma)
			{
				cw[0]=0.0;
				ccw[0]=p[1].gamma;
				cw[1]=p[0].gamma;
				ccw[1]=PIt2;
				nnr=2;
			}
			else
			{
				cw[0]=p[0].gamma;
				ccw[0]=p[1].gamma;
				nnr=1;
			}
		}
		else 
		{
			if(p[0].gamma>p[1].gamma)
			{
				cw[0]=0.0;
				ccw[0]=p[1].gamma;
				cw[1]=p[0].gamma;
				ccw[1]=PIt2;
				cw[2]=p[2].gamma;
				ccw[2]=p[3].gamma;
				nnr=3;
			}
			else
			{
				cw[0]=p[0].gamma;
				ccw[0]=p[1].gamma;
				if(p[2].gamma>p[3].gamma)
				{
					cw[1]=p[2].gamma;
					ccw[1]=PIt2;
					cw[2]=0.0;
					ccw[2]=p[3].gamma;
					nnr=3;
				}
				else
				{
					cw[1]=p[2].gamma;
					ccw[1]=p[3].gamma;
					nnr=2;
				}
			}
		}		//region made

		for(i=0;i<nnr;i++)
		{
			cwr=(-1);
			ccwr=(-1);
			for(j=0;j<enreg[sn];j++)	//what existing region are the end points new of region i in?
			{
				if(cw[i]>=egamma[sn][2*j]&&cw[i]<=egamma[sn][2*j+1])
					cwr=j;
				if(ccw[i]>=egamma[sn][2*j]&&ccw[i]<=egamma[sn][2*j+1])
					ccwr=j;
			}
			if(cwr<0&&ccwr<0)	//both are not in an existing region 
			{
				k=(-1);
				for(j=0;j<enreg[sn]&&k<0;j++)
					if(egamma[sn][2*j]>=cw[i]&&egamma[sn][2*j]<=ccw[i])
						k=j;
				if(k>=0)	//new region contains an existing region
				{
					egamma[sn][2*k]=cw[i];
					egamma[sn][2*k+1]=ccw[i];
				}
				else   //or it doesn't
				{
					if(enreg[sn]>=3)
						return -1;
					egamma[sn][2*enreg[sn]]=cw[i];
					egamma[sn][2*enreg[sn]+1]=ccw[i];
					enreg[sn]++;
				}
			}
			else if(cwr>=0&&ccwr>=0&&cwr!=ccwr)		//they are in different regions
			{
				egamma[sn][2*cwr+1]=egamma[sn][2*ccwr+1];
				egamma[sn][2*ccwr]=egamma[sn][2*cwr];
				if(enreg[sn]==3)
					if((cwr==0&&ccwr==1)||(cwr==1&&ccwr==0))
					{
						egamma[sn][2]=egamma[sn][4];
						egamma[sn][3]=egamma[sn][5];
					}
				enreg[sn]--;
			}
			else if(cwr<0||ccwr<0)			//one is in a region, one is not
			{
				if(cwr>=0)
					egamma[sn][2*cwr+1]=ccw[i];
				else
					egamma[sn][2*ccwr]=cw[i];
			}
			cwr=(-1);		//check if one region contains another
			ccwr=(-1);
			for(j=0;j<enreg[sn]-1;j++)
				for(k=j+1;k<enreg[sn];k++)
				{
					if(egamma[sn][2*j]>egamma[sn][2*k]&&egamma[sn][2*j]<egamma[sn][2*k+1])	//k contains j
						if(egamma[sn][2*j+1]>egamma[sn][2*k]&&egamma[sn][2*j+1]<egamma[sn][2*k+1])
						{
							ccwr=k;
							cwr=j;
						}
						else
							return -2;
					if(egamma[sn][2*k]>egamma[sn][2*j]&&egamma[sn][2*k]<egamma[sn][2*j+1])	//j contains k
						if(egamma[sn][2*k+1]>egamma[sn][2*j]&&egamma[sn][2*k+1]<egamma[sn][2*j+1])
						{
							ccwr=j;
							cwr=k;
						}
						else return -2;
					if(cwr>=0)
					{
						if(enreg[sn]==2)
						{
							if(ccwr==1&&cwr==0)
							{
								egamma[sn][0]=egamma[sn][2];
								egamma[sn][1]=egamma[sn][3];
							}
							enreg[sn]--;
						}
						else if(enreg[sn]==3)
						{
							if(cwr!=2)
							{
								egamma[sn][2*cwr]=egamma[sn][4];
								egamma[sn][2*cwr+1]=egamma[sn][5];
							}
							enreg[sn]--;
						}
						else
							return -3;
						j=(-1);
						k=0;
					}
				}
		}
	}
	else if(eorc=='c')
	{
		if(np=2&&ceregactive[sn][1]>=0)
		{
			cw[0]=p[0].gamma;
			ccw[0]=p[1].gamma;
			if(cw[0]>ccw[0])
				return -4;
			if(cw[0]>0.0&&ccw[0]<PI)
			{
				if(ceregactive[sn][1])
				{
					if(cw[0]<cegamma[sn][1])
						cegamma[sn][1]=cw[0];
					if(ccw[0]>cegamma[sn][2])
						cegamma[sn][2]=ccw[0];
				}
				else
				{
					cegamma[sn][1]=cw[0];
					cegamma[sn][2]=ccw[0];
					ceregactive[sn][1]=1;
				}
			}
			else if(cw[0]==0.0)
			{
				if(ceregactive[sn][0])
				{
					if(ccw[0]>cegamma[sn][0])
						cegamma[sn][0]=ccw[0];
				}
				else
				{
					cegamma[sn][0]=ccw[0];
					ceregactive[sn][0]=1;
				}
			}
			else if(ccw[0]==PI)
			{
				if(ceregactive[sn][2])
				{
					if(cw[0]<cegamma[sn][3])
						cegamma[sn][3]=cw[0];
				}
				else
				{
					cegamma[sn][3]=cw[0];
					ceregactive[sn][2]=1;
				}
			}
			if(ceregactive[sn][1]&&ceregactive[sn][0])
				if(cegamma[sn][0]>cegamma[sn][1])
				{
					cegamma[sn][0]=cegamma[sn][2];
					ceregactive[sn][1]=0;
				}
			if(ceregactive[sn][1]&&ceregactive[sn][2])
				if(cegamma[sn][3]<cegamma[sn][2])
				{
					cegamma[sn][3]=cegamma[sn][1];
					ceregactive[sn][1]=0;
				}
			if(ceregactive[sn][0]&&ceregactive[sn][2])
				if(cegamma[sn][0]>ceregactive[sn][3])
					ceregactive[sn][1]=(-1);	//signify that whole arc is covered
		}
		if(nq==2)
		{
			cw[0]=q[0].gamma;
			ccw[0]=q[1].gamma;
			if(cw[0]>ccw[0])
				return -5;
			if(cw[0]>0.0&&ccw[0]<PI)
			{
				if(ccregactive[sn][1])
				{
					if(cw[0]<ccgamma[sn][1])
						ccgamma[sn][1]=cw[0];
					if(ccw[0]>ccgamma[sn][2])
						ccgamma[sn][2]=ccw[0];
				}
				else
				{
					ccgamma[sn][1]=cw[0];
					ccgamma[sn][2]=ccw[0];
					ccregactive[sn][1]=1;
				}
			}
			else if(cw[0]==0.0)
			{
				if(ccregactive[sn][0])
				{
					if(ccw[0]>ccgamma[sn][0])
						ccgamma[sn][0]=ccw[0];
				}
				else
				{
					ccgamma[sn][0]=ccw[0];
					ccregactive[sn][0]=1;
				}
			}
			else if(ccw[0]==PI)
			{
				if(ccregactive[sn][2])
				{
					if(cw[0]<ccgamma[sn][3])
						ccgamma[sn][3]=cw[0];
				}
				else
				{
					ccgamma[sn][3]=cw[0];
					ccregactive[sn][2]=1;
				}
			}
			if(ccregactive[sn][1]&&ccregactive[sn][0])
				if(ccgamma[sn][0]>ccgamma[sn][1])
				{
					ccgamma[sn][0]=ccgamma[sn][2];
					ccregactive[sn][1]=0;
				}
			if(ccregactive[sn][1]&&ccregactive[sn][2])
				if(ccgamma[sn][3]<ccgamma[sn][2])
				{
					ccgamma[sn][3]=ccgamma[sn][1];
					ccregactive[sn][1]=0;
				}
			if(ccregactive[sn][0]&&ccregactive[sn][2])
				if(ccgamma[sn][0]>ccregactive[sn][3])
					ccregactive[sn][1]=(-1);	//signify that whole arc is covered
		}
	}
	return 0;
}
int areacovered(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,char acoutfilename[128])
{                       //finds portion of each spot that is never covered by the planet during one pass -- not perfect, but a reasonable approximation when it works
	int i,j,k,np,nq,wp[4],wq[2],en,ti,a,b;
	double starradius,t,ydif,zdif,dsq,dmax,vc,wc,lambda,y,z,dang;
	vwpoint p[4],q[4];
	double egamma[MAXSPOTS][6],ccgamma[MAXSPOTS][4],cegamma[MAXSPOTS][4];
	int ebeencovered[MAXSPOTS],cbeencovered[MAXSPOTS],enreg[MAXSPOTS],ccregactive[MAXSPOTS][3],ceregactive[MAXSPOTS][3];
	double thisarea[MAXSPOTS],maxarea[MAXSPOTS];
	int tindexofmaxarea[MAXSPOTS];
	double eareauncovered,careauncovered,totalarea,caparea;
	double scaledcoveredarea[MAXSPOTS];

	FILE *acout;

	starradius=star->r;
	errorflag=0;

	for(i=0;i<numspots;i++)
	{
		ebeencovered[i]=0; cbeencovered[i]=0;
		enreg[i]=0; 
		maxarea[i]=0;
		for(j=0;j<3;j++)
		{
			ccregactive[i][j]=0;
			ceregactive[i][j]=0;
		}
	}

	star->phi=star->omega*lctime[(int)(lcn/2)];

	acout=fopen(acoutfilename,"w");
	for(ti=0;ti<lcn;ti++)
	{
		t=lctime[ti];
		zerototalnums();

		for(i=0;i<numplanets;i++)
		{
			setplanetpos(t,planet,i);
			if(planet[i].x>0)	//planet is in front of star
				if(sqrt(planet[i].y*planet[i].y+planet[i].z*planet[i].z)-planet[i].r<star->r)
					createcircle(planet[i],star->r,star->rsq);
		}

		updatespots(spot,star->r,star->phi);
		for(i=0;i<numspots;i++)
			if(spot[i].psi-spot[i].alpha<PIo2)	//spot visable
				createellipsecresent(i,spot[i],star->r);

		if(totalnum.circle>0)
		{
			for(i=0;i<totalnum.ellipse;i++)	//-----------check ellipses
			{
				ellipse[i].numcircleint[0]=0;
				ydif=circle[0].centery-ellipse[i].centery;
				zdif=circle[0].centerz-ellipse[i].centerz;
				dsq=ydif*ydif+zdif*zdif;
				dmax=circle[0].radius+ellipse[i].smajor;
				if(dsq<dmax*dmax)
				{
					vc=ydif*ellipse[i].vhaty+zdif*ellipse[i].vhatz;	//v-w coords of planet center
					wc=ydif*ellipse[i].whaty+zdif*ellipse[i].whatz;	//yes, really
					np=ellipsecircleintersection(ellipse[i].smajor,ellipse[i].sminor,vc,wc,circle[0].rsq,p);
					if(np==0)
					{
						if(ellipse[i].active)
						{
							if((vc-ellipse[i].smajor)*(vc-ellipse[i].smajor)+wc*wc<circle[0].rsq)	//point on ellipse inside circle
							{
								ellipse[i].active=0;	//ellipse completely covered
								ebeencovered[ellipse[i].spot]=1;
							}
							else if((vc+circle[0].radius)*(vc+circle[0].radius)/(ellipse[i].smajor*ellipse[i].smajor)+wc*wc/(ellipse[i].sminor*ellipse[i].sminor)<1.0)	//point on circle inside ellipse
							{	
								ellipsehole[totalnum.ellipsehole].ellipsen=i;
								ellipsehole[totalnum.ellipsehole].circlen=0;
								ellipsehole[totalnum.ellipsehole].area=ellipse[i].area-circle[0].area;
								totalnum.ellipsehole++;
								ellipse[i].active=0;
							}
						}
					}
					else if(np==2)
					{
						ellipse[i].numcircleint[0]=2;
						for(j=0;j<2;j++)
						{
							ellipse[i].circleint[0][j].v=p[j].v;	//ellipse under circle from 0 to 1
							ellipse[i].circleint[0][j].w=p[j].w;
							ellipse[i].circleint[0][j].gamma=p[j].gamma;	
						}
						if(ellipse[i].active)
						{
							createellipsepart(p[1],p[0],i,0,vc,wc);
							ellipse[i].active=0;
							errorflag=areacoveredupdate(ellipse[i].spot,'e',2,0,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
						}
					}
					else if(np==4)
					{
						ellipse[i].numcircleint[0]=4;
						for(j=0;j<4;j++)
						{
							ellipse[i].circleint[0][j].v=p[j].v;	//ellipse under circle: 0 to 1, 2 to 3
							ellipse[i].circleint[0][j].w=p[j].w;
							ellipse[i].circleint[0][j].gamma=p[j].gamma;
						}
						if(ellipse[i].active)
						{
							createellipsepart(p[1],p[2],i,0,vc,wc);
							createellipsepart(p[3],p[0],i,0,vc,wc);
							ellipse[i].active=0;
							errorflag=areacoveredupdate(ellipse[i].spot,'e',4,0,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
						}
					}
					else
						errorflag=4;
				}
			}


			for(i=0;i<totalnum.cresent;i++)	////-----------check cresents
			{
				k=0;
				en=cresent[i].ellipsen;
				np=0;	//find circle intersections with ellipse arc
				for(j=0;j<ellipse[en].numcircleint[0];j++)
					if(isanglebetween(ellipsearc[cresent[i].ellipsearcn].end[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma,ellipse[en].circleint[0][j].gamma))
					{
						wp[np]=j;
						np++;
					}
				nq=0;	//find circle intersections with scircle arc
				for(j=0;j<circle[0].numstarint;j++)
					if(isanglebetween(scirclearc[cresent[i].scirclearcn].ebeta[0],scirclearc[cresent[i].scirclearcn].ebeta[1],circle[0].starintbeta[j]))
					{
						if(nq>=2)
							errorflag=4;
						wq[nq]=j;
						nq++;
					}
				if(np%2!=nq%2||np>2)
				{
					if(!fixcresentnpnq(starradius,i,&np,&nq,wp,wq))
						errorflag=4;
				}

				ydif=circle[0].centery-ellipse[en].centery;
				zdif=circle[0].centerz-ellipse[en].centerz;
				vc=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	//v-w coords of planet center
				wc=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	//yes, really
				if(nq==0)
				{
					if(np==0)
					{
						if((ellipsearc[cresent[i].ellipsearcn].end[0].v-vc)*(ellipsearc[cresent[i].ellipsearcn].end[0].v-vc)+(ellipsearc[cresent[i].ellipsearcn].end[0].w-wc)*(ellipsearc[cresent[i].ellipsearcn].end[0].w-wc)<circle[0].rsq)
						{
							cresent[i].active=0;	//cresent covered by circle
							cbeencovered[ellipse[cresent[i].ellipsen].spot]=1;
						}
						else if(circle[0].mindcen>cresent[i].mindcen&&circle[0].maxdcen<scirclearc[cresent[i].scirclearcn].radius
								&&isanglebetween(ellipsearc[cresent[i].ellipsearcn].end[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma,atan2(wc,vc)))
						{	//circle is completly inside cresent
							cresenthole[totalnum.cresenthole].circlen=0;
							cresenthole[totalnum.cresenthole].cresentn=i;
							cresenthole[totalnum.cresenthole].area=cresent[i].area-circle[0].area;
							totalnum.cresenthole++;
							cresent[i].active=0;
							whichspots[ellipse[cresent[i].ellipsen].spot]=1;
						}
					}
					else if(np==2)	//np=2 nq=0
					{
						createcresentpart20(ellipse[en].circleint[0][wp[0]],ellipse[en].circleint[0][wp[1]],i,en,0,vc,wc);
						cresent[i].active=0;
						errorflag=areacoveredupdate(ellipse[cresent[i].ellipsen].spot,'c',2,0,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
					}
					else
						errorflag=4;
				}
				else if(nq==2)
				{
					for(j=0;j<2;j++)
					{
						ydif=starradius*cos(circle[0].starintbeta[wq[j]])-ellipse[en].centery;
						zdif=starradius*sin(circle[0].starintbeta[wq[j]])-ellipse[en].centerz;
						q[j].v=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	
						q[j].w=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	
						q[j].gamma=atan2(q[j].w,q[j].v);
						if(q[j].gamma<0) q[j].gamma+=PIt2;
					}
					if(np==0)	//np=0 nq=2
					{
						createcresentpart02(q[0],q[1],i,en,0,vc,wc,starradius);
						cresent[i].active=0;
						errorflag=areacoveredupdate(ellipse[cresent[i].ellipsen].spot,'c',0,2,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
					}
					else if(np==2)	//np=2 nq=2
					{
						createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],1,i,en,0,vc,wc,starradius);
						createcresentpart11(ellipse[en].circleint[0][wp[1]],q[1],0,i,en,0,vc,wc,starradius);
						cresent[i].active=0;
						errorflag=areacoveredupdate(ellipse[cresent[i].ellipsen].spot,'c',2,2,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
					}
					else
						errorflag=4;
				}
				else if(nq==1)
				{
					ydif=starradius*cos(circle[0].starintbeta[wq[0]])-ellipse[en].centery;
					zdif=starradius*sin(circle[0].starintbeta[wq[0]])-ellipse[en].centerz;
					q[0].v=ydif*ellipse[en].vhaty+zdif*ellipse[en].vhatz;	
					q[0].w=ydif*ellipse[en].whaty+zdif*ellipse[en].whatz;	
					q[0].gamma=atan2(q[0].w,q[0].v);
					if(q[0].gamma<0) q[0].gamma+=PIt2;
					if(np==1)
					{
						if(wq[0]%2!=wp[0]%2)
							errorflag=4;
						if(wq[0]==0)
						{                 //clockwise side
							createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],1,i,en,0,vc,wc,starradius);
							p[0].gamma=0.0;  //not really zero, just the end
							p[1].gamma=ellipse[en].circleint[0][wp[0]].gamma;
							q[1].gamma=q[0].gamma;
							q[0].gamma=0.0;
						}
						else
						{                 //counterclockwise side
							createcresentpart11(ellipse[en].circleint[0][wp[0]],q[0],0,i,en,0,vc,wc,starradius);
							p[0].gamma=ellipse[en].circleint[0][wp[0]].gamma;
							p[1].gamma=PI;  //not really PI, just the end
							q[1].gamma=PI;
						}
						errorflag=areacoveredupdate(ellipse[cresent[i].ellipsen].spot,'c',2,2,p,q,egamma,ccgamma,cegamma,ebeencovered,cbeencovered,enreg,ccregactive,ceregactive);
						cresent[i].active=0;
					}
					else
						errorflag=4;
				}
				else
					errorflag=4;
			}
		}
		if(errorflag)
			return errorflag;
		for(i=0;i<numspots;i++)
			thisarea[i]=0;
		for(i=0;i<totalnum.ellipse;i++)
			thisarea[ellipse[i].spot]=ellipse[i].area;
		for(i=0;i<totalnum.cresent;i++)
			thisarea[ellipse[cresent[i].ellipsen].spot]+=cresent[i].area;
		for(i=0;i<numspots;i++)
			if(thisarea[i]>maxarea[i])
			{
				maxarea[i]=thisarea[i];
				tindexofmaxarea[i]=ti;
			}
	}

	for(k=0;k<numspots;k++)
	{
		t=lctime[tindexofmaxarea[k]];
		zerototalnums();

		star->phi=star->omega*t;

		for(i=0;i<numplanets;i++)
		{
			setplanetpos(t,planet,i);
			if(planet[i].x>0)	//planet is in front of star
				if(sqrt(planet[i].y*planet[i].y+planet[i].z*planet[i].z)-planet[i].r<star->r)
					createcircle(planet[i],star->r,star->rsq);
		}

		updatespots(spot,star->r,star->phi);
		for(i=0;i<numspots;i++)
			if(spot[i].psi-spot[i].alpha<PIo2)	//spot visable
				createellipsecresent(i,spot[i],star->r);

		totalarea=0.0;
		eareauncovered=0.0;
		i=(-1);
		for(j=0;j<totalnum.ellipse;j++)
			if(ellipse[j].spot==k)
				i=j;
		if(i>=0)	//ellipse i is spot k
		{
			totalarea=ellipse[i].area;
			if(ebeencovered[k])		//totally covered
				eareauncovered=0.0;
			else if(enreg[k]==0)	//totally uncovered
				eareauncovered=ellipse[i].area;
			else					//somewhat covered
			{
				if(enreg[k]>1)	//if split at zero, combine
				{
					a=(-1);
					b=(-1);
					for(j=0;j<enreg[k];j++)
					{
						if(egamma[k][2*j]==0.0)
							b=j;
						if(egamma[k][2*j+1]==PIt2)
							a=j;
					}
					if(a>=0&&b>=0)
					{
						egamma[k][2*b]=egamma[k][2*a];
						egamma[k][2*a+1]=egamma[k][2*b+1];
						if(a<2&&b<2)
						{
							egamma[k][2]=egamma[k][4];
							egamma[k][3]=egamma[k][5];
						}
						enreg[k]--;
					}
				}
				
				if(enreg[k]==1)
					eareauncovered=areaellipsesection(ellipse[i].smajor,ellipse[i].sminor,egamma[k][1],egamma[k][0]);
				else if(enreg[k]==2)
					eareauncovered=areaellipsesection(ellipse[i].smajor,ellipse[i].sminor,egamma[k][1],egamma[k][2])+areaellipsesection(ellipse[i].smajor,ellipse[i].sminor,egamma[k][3],egamma[k][0]);
				else
					return 5;
			}
		}		//end of ellipse part

		careauncovered=0.0;
		i=(-1);
		for(j=0;j<totalnum.cresent;j++)
			if(ellipse[cresent[j].ellipsen].spot==k)
				i=j;
		if(i>=0)	//cresent i is spot k
		{
			totalarea+=cresent[i].area;
			if(cbeencovered[k])
				careauncovered=0.0;
			else
			{
				a=0; b=0;
				for(j=0;j<3;j++)
				{
					if(ccregactive[k][j])
						a++;
					if(ceregactive[k][j])
						b++;
				}

				if(a>2||b>2)
					return 6;

				if(a==0)
				{
					if(b!=0)
						return 7;
					if(ccregactive[k][1]<0&&ceregactive[k][1]<0)		//both sides all covered
						careauncovered=0.0;	
					else if(ccregactive[k][1]==0&&ceregactive[k][1]==0)	//both sides all uncovered
						careauncovered=cresent[i].area;
					else
						return 8;
				}
				else
				{
					p[0].gamma=p[1].gamma=p[2].gamma=q[0].gamma=q[1].gamma=q[2].gamma=(-1.0);	//set up p and q 
					if(ceregactive[k][0])
					{
						p[0].gamma=cegamma[k][0];
						if(ceregactive[k][1])
						{
							p[1].gamma=cegamma[k][1];
							p[2].gamma=cegamma[k][2];
							if(ceregactive[k][2])
								return 9;
						}
						else if(ceregactive[k][2])
							p[1].gamma=cegamma[k][3];
					}
					else
					{
						if(ceregactive[k][1])
						{
							p[0].gamma=cegamma[k][1];
							p[1].gamma=cegamma[k][2];
							if(ceregactive[k][2])
								p[2].gamma=cegamma[k][3];
						}
						else if(ceregactive[k][2])
							p[0].gamma=cegamma[k][3];
					}
					if(ccregactive[k][0])
					{
						q[0].gamma=ccgamma[k][0];
						if(ccregactive[k][1])
						{
							q[1].gamma=ccgamma[k][1];
							q[2].gamma=ccgamma[k][2];
							if(ccregactive[k][2])
								return 10;
						}
						else if(ccregactive[k][2])
							q[1].gamma=ccgamma[k][3];
					}
					else
					{
						if(ccregactive[k][1])
						{
							q[0].gamma=ccgamma[k][1];
							q[1].gamma=ccgamma[k][2];
							if(ccregactive[k][2])
								q[2].gamma=ccgamma[k][3];
						}
						else if(ccregactive[k][2])
							q[0].gamma=ccgamma[k][3];
					}
					for(j=0;j<3;j++)
					{
						if(p[j].gamma>=0)
						{
							lambda=atan(ellipse[cresent[i].ellipsen].smajor*tan(p[j].gamma)/ellipse[cresent[i].ellipsen].sminor);
							if(p[j].gamma>PIo2&&p[j].gamma<thPIo2)
								lambda+=PI;
							p[j].v=ellipse[k].smajor*cos(lambda);
							p[j].w=ellipse[k].sminor*sin(lambda);
						}
						if(q[j].gamma>=0)
						{
							lambda=atan(ellipse[cresent[i].ellipsen].smajor*tan(q[j].gamma)/ellipse[cresent[i].ellipsen].sminor);
							if(q[j].gamma>PIo2&&q[j].gamma<thPIo2)
								lambda+=PI;
							q[j].v=ellipse[k].smajor*cos(lambda);
							q[j].w=ellipse[k].sminor*sin(lambda);
						}
					}		//p and q are set

					if(b==0)
					{
						if(ccregactive[k][0]==ccregactive[k][2]&&ccregactive[k][0]!=ccregactive[k][1])
						{
							vwtoyz(p[1].v,p[1].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							vwtoyz(p[0].v,p[0].w,cresent[i].ellipsen,&y,&z);
							dang=dang-atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered=areacirclesection(starradius,dang);
							if(ccregactive[k][1])
								careauncovered=cresent[i].area-careauncovered;
						}
						else
							return 9;
					}
					else if(a==1&&b==1&&!ccregactive[k][1])
					{
						if(ccregactive[k][0]&&ceregactive[k][0])
						{
							careauncovered=areatriangle(p[0],q[0],ellipsearc[cresent[i].ellipsearcn].end[1]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
							dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
						}
						else if(ccregactive[k][2]&&ceregactive[k][2])
						{
							careauncovered=areatriangle(q[0],p[0],ellipsearc[cresent[i].ellipsearcn].end[0]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[0].gamma);
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
						}
						else
							return 10;
					}
					else if((a==1&&b==1&&ccregactive[k][1]&&ceregactive[k][1])||(a==2&&b==2&&(!ccregactive[k][1])&&(!ceregactive[k][1])))
					{
						careauncovered=areatriangle(p[1],q[1],ellipsearc[cresent[i].ellipsearcn].end[1]);
						careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[1].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
						dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
						vwtoyz(q[1].v,q[1].w,cresent[i].ellipsen,&y,&z);
						dang-=atan2(z,y);
						if(dang<0)
							dang+=PIt2;
						careauncovered+=areacirclesection(starradius,dang);

						careauncovered+=areatriangle(q[0],p[0],ellipsearc[cresent[i].ellipsearcn].end[0]);
						careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[0].gamma);
						vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
						dang=atan2(z,y);
						dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
						if(dang<0)
							dang+=PIt2;
						careauncovered+=areacirclesection(starradius,dang);

						if(a==2)
							careauncovered=cresent[i].area-careauncovered;
					}
					else if(a==2&&b==1&&(!ceregactive[k][1]))
					{
						if(ccregactive[k][0]&&ceregactive[k][0])
						{
							careauncovered=areatriangle(p[0],q[2],ellipsearc[cresent[i].ellipsearcn].end[1]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[0].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
							dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
							vwtoyz(q[2].v,q[2].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
							vwtoyz(q[1].v,q[1].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
						}
						else if(ccregactive[k][2]&&ceregactive[k][2])
						{
							careauncovered=areatriangle(q[0],p[0],ellipsearc[cresent[i].ellipsearcn].end[0]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[0].gamma);
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
							vwtoyz(q[2].v,q[2].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							vwtoyz(q[1].v,q[1].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);
						}
						else
							return 11;

					}
					else if(a==2&&b==2&&ccregactive[k][1]&&ceregactive[k][1])
					{
						if((!ceregactive[k][2])&&(!ccregactive[k][2]))
						{
							careauncovered=areatriangle(p[2],q[2],ellipsearc[cresent[i].ellipsearcn].end[1]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[2].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
							dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
							vwtoyz(q[2].v,q[2].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);

							careauncovered+=areatriangle(q[1],p[1],ellipsearc[cresent[i].ellipsearcn].end[0]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[1].gamma);
							vwtoyz(q[1].v,q[1].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);

							careauncovered-=areatriangle(q[0],p[0],ellipsearc[cresent[i].ellipsearcn].end[0]);
							careauncovered+=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[0].gamma);
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
							if(dang<0)
								dang+=PIt2;
							careauncovered-=areacirclesection(starradius,dang);
						}
						else if((!ceregactive[k][0])&&(!ccregactive[k][0]))
						{
							careauncovered=areatriangle(p[1],q[1],ellipsearc[cresent[i].ellipsearcn].end[1]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[1].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
							dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
							vwtoyz(q[1].v,q[1].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);

							careauncovered+=areatriangle(q[0],p[0],ellipsearc[cresent[i].ellipsearcn].end[0]);
							careauncovered-=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,ellipsearc[cresent[i].ellipsearcn].end[0].gamma,p[0].gamma);
							vwtoyz(q[0].v,q[0].w,cresent[i].ellipsen,&y,&z);
							dang=atan2(z,y);
							dang-=scirclearc[cresent[i].scirclearcn].ebeta[0];
							if(dang<0)
								dang+=PIt2;
							careauncovered+=areacirclesection(starradius,dang);

							careauncovered-=areatriangle(p[2],q[2],ellipsearc[cresent[i].ellipsearcn].end[1]);
							careauncovered+=areaellipsesection(ellipse[cresent[i].ellipsen].smajor,ellipse[cresent[i].ellipsen].sminor,p[2].gamma,ellipsearc[cresent[i].ellipsearcn].end[1].gamma);
							dang=scirclearc[cresent[i].scirclearcn].ebeta[1];
							vwtoyz(q[2].v,q[2].w,cresent[i].ellipsen,&y,&z);
							dang-=atan2(z,y);
							if(dang<0)
								dang+=PIt2;
							careauncovered-=areacirclesection(starradius,dang);
						}
						else
							return 12;
					}
					else  
						return 13;
				}
			}
		}	//end of cresent part

		fprintf(acout,"spot %i\n projected area covered:%lf\n total projected area: %lf\n",k,totalarea-(eareauncovered+careauncovered),totalarea);
		printf("spot %i\n projected area covered:%lf\n total projected area: %lf\n",k,totalarea-(eareauncovered+careauncovered),totalarea);

		caparea=PIt2*starradius*starradius*(1.0-cos(spot[k].alpha));
		if(totalarea>0.0)
			scaledcoveredarea[k]=(totalarea-eareauncovered-careauncovered)/totalarea;	//fraction covered 
		else
			scaledcoveredarea[k]=0.0;

		fprintf(acout,"     projected fraction covered: %lf\n total area on sperical surface: %lf\n",scaledcoveredarea[k],caparea);
		printf("     projected fraction covered: %lf\n total area on sperical surface: %lf\n",scaledcoveredarea[k],caparea);

		scaledcoveredarea[k]*=caparea;	//times area of spherical cap

		fprintf(acout,"         total covered area on surface of sphere: %lf\n\n",scaledcoveredarea[k]);
		printf("         total covered area on surface of sphere: %lf\n\n",scaledcoveredarea[k]);
	}

	caparea=PI*starradius*starradius*(2.0*planet[0].r);	//area of ribbon covered by planet
	totalarea=0.0;
	for(i=0;i<numspots;i++)
		totalarea+=scaledcoveredarea[i];

	fprintf(acout,"\n       sum of spot areas covered: %lf\narea of ribbon covered by planet: %lf\n\nfraction of ribbon with spots: %lf\n",totalarea,caparea,totalarea/caparea);
	printf("\n       sum of spot areas covered: %lf\narea of ribbon covered by planet: %lf\n\nfraction of ribbon with spots: %lf\n",totalarea,caparea,totalarea/caparea);
	
	fclose(acout);
	return 0;
}
void debuglc(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[])
{
	int i;
	double light,*param;
	FILE *lcout;
	
	PRINTVIS=WHICHPRINTVIS;
	pvasc=1; pvisc=10;
#	if ANYPRINTVIS
		outv=fopen("vis.txt","w");
#	endif

	param=(double *)malloc(3*numspots*sizeof(double));
		if(param==NULL)
			printf("malloc error\n");

	for(debugtrial=0;debugtrial<4;debugtrial++)
		setrandomparam(param,star);
	for(;debugtrial<5000;debugtrial++)
	{
		printf("trial: %i\n",debugtrial);
		lcout=fopen("lcout.txt","w");
		setrandomparam(param,star);
		(void)setspots(param,spot,star,planet);
		for(i=0;i<lcn;i++)
		{
			debuglcni=i;
			errorflag=0;
			light=lightness(lctime[i],star,planet,spot);
			if(!(light>=0&&light<=3.1416))
				printf("light error\n");
			fprintf(lcout,"%lf %lf\n",lclight[i],light);
			if(errorflag)
			{
				fprintf(outerr,"error %i at datum %i (location 2)\n",errorflag,i);
				printf("debuglc: error\n");
			}
#			if ANYPRINTVIS
				if(PRINTVIS)
				{
					fprintf(outv,"z\n%lf %lf\n",lclight[i],light);
					fflush(outv);
				}
#			endif
		}
		fclose(lcout);
	}
	free((void *)param);
}
double degenfindr(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],double th,double light0)
{
	double rmin=0.0,rmax=.99,light,dl=100.0,param[3];

	param[1]=th;
	param[2]=0;
	while(dl>0.0000001)
	{
		param[0]=(rmin+rmax)/2.0;
		setspots(param,spot,star,planet);
		light=lightness(0.0,star,planet,spot);
		if(light>light0)
		{
			rmin=param[0];
			dl=light-light0;
		}
		else
		{
			rmax=param[0];
			dl=light0-light;
		}
	}	
	return param[0];
}
void degen(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS])
{	//just for investigating latitude-size degeneracy 
	char filename[64];
	int i,noc=20;
	double r0,dth,param[3],light0;
	double t,dt,halfperiod;
	FILE *out;

	r0=0.1;
	param[0]=r0;
	param[1]=PIo2;
	param[2]=0.0;
	setspots(param,spot,star,planet);
	dth=(PIo2)/((double)(noc-1));
	light0=lightness(0.0,star,planet,spot);
	halfperiod=5.0*24.0*60.0*60.0;
	dt=halfperiod/5000.0;

	for(i=0;i<noc;i++)
	{
		sprintf(filename,"degenout%i.txt",i);
		out=fopen(filename,"w");
		param[1]=PIo2-dth*i;
		param[0]=degenfindr(star,planet,spot,param[1],light0);
		fprintf(out,"%0.12lf\n%0.12lf\n",param[0],param[1]);
		setspots(param,spot,star,planet);
		for(t=(-halfperiod);t<=halfperiod;t+=dt)
			fprintf(out,"%0.12lf\n",lightness(t,star,planet,spot));
		fclose(out);
	}
}
void gentimetest(stardata *star,planetdata planet[MAXPLANETS],spotdata spot[MAXSPOTS],int lcn,double lctime[],double lclight[],double lcuncertainty[],double lclightnorm,char tocheck,double readparam[])
{		//produce model lightcurve with time varying spots
	int i,j;
	double torig,light;
	double param[15],rmax,t0[5],tmax[5],tf[5],fr[5],fd[5];
	double decayf=0.0000000005; //0.0000000002588;	//decay factor
	double risetf=10.0;	//rise time factor (rise time = 1/risetf of decay time)
	double decayfrat=100.0; //ratio involved with decay time
	double dfoPI;
	FILE *lcout;

	dfoPI=decayfrat/PI;
	PRINTVIS=WHICHPRINTVIS;
	pvasc=2; pvisc=4;
#	if ANYPRINTVIS
		outv=fopen("stimetest_vis.txt","a");
#	endif

	setrandomparam(param,star);
	for(i=0;i<5;i++)
	{
		rmax=param[3*i];
		t0[i]=lctime[0]-5000.0+RND*(lctime[lcn-1]-lctime[0])/2.0;
		tmax[i]=t0[i]-log(decayfrat/(PI*rmax*rmax+decayfrat))/(risetf*decayf);
		tf[i]=tmax[i]-log(decayfrat/(PI*rmax*rmax+decayfrat))/decayf;
		fr[i]=rmax/sqrt(tmax[i]-t0[i]);
		fd[i]=rmax*rmax+dfoPI;
	}


	if(!tocheck)
	{
		lcout=fopen("timetestgoal.txt","w");
		fprintf(lcout,"overall time %lf  -  %lf\n",lctime[0]/86400.0+planet[0].lct0,lctime[lcn-1]/86400.0+planet[0].lct0);
		for(i=0;i<5;i++)
		{
			fprintf(lcout,"\nspot %i:\nrmax= %lf\ntheta= %lf\n phi=%lf\n",i,param[i*3],param[i*3+1],param[i*3+2]);
			fprintf(lcout,"time %lf  %lf  %lf\n",t0[i]/86400.0+planet[0].lct0,tmax[i]/86400.0+planet[0].lct0,tf[i]/86400.0+planet[0].lct0);
		}
		fclose(lcout);
	
		lcout=fopen("timetest.lc","w");
	}

	for(i=0;i<lcn;i++)
	{
		debuglcni=i;
		errorflag=0;
		torig=lctime[i]/86400.0+planet[0].lct0;
		for(j=0;j<numspots;j++)
		{
			if(lctime[i]<t0[j])
				param[j*3]=0.0;
			else if(lctime[i]<tmax[j])
				param[j*3]=fr[j]*sqrt(lctime[i]-t0[j]);
			else if(lctime[i]<tf[j])
				param[j*3]=sqrt(fd[j]*exp(decayf*(tmax[j]-lctime[i]))-dfoPI);
			else
				param[j*3]=0.0;
			if(!(param[j*3]>TOOSMALL&&param[j*3]<TOOBIG))
				printf("pig\n");
		}
		setspots(param,spot,star,planet);
		light=lightness(lctime[i],star,planet,spot);
		if(!tocheck)
			fprintf(lcout,"%0.9lf %0.6lf %0.6lf\n",torig,light/lclightnorm,lcuncertainty[i]/lclightnorm);	//for making test lightcurves
		if(errorflag)
		{
			fprintf(outerr,"error %i at datum %i (location 3)\n",errorflag,i);
			printf("lcgen: error\n");
			exit(0);
		}
		if(tocheck)
		{
			lclight[i]=light;
			setspots(readparam,spot,star,planet);
			pvasc=1; pvisc=10;
			light=lightness(lctime[i],star,planet,spot);
			pvasc=2; pvisc=4;
		}
#		if ANYPRINTVIS
			if(PRINTVIS)
			{
				if(PRINTVIS!=1)
					fprintf(outv,"z\n%lf %lf\n",lclight[i],light);
			}
#		endif
	}
	if(!tocheck)
		fclose(lcout);

}

int main(int argc,char *argv[])
{
	char filename[64],rootname[64],seedfilename[64];
	int i,j;
	int lcn;
	int mcmcnpop,randomseed;
	int partitionpop,partitionsteps;
	long int mcmcmaxstepsortime;
	long int memused[2];
	double *lctime,*lclight,*lcuncertainty,lclightnorm,lcstarttime,lcfinishtime,lcmaxlight,ascale;
	double *readparam;
	stardata *star,thestar;
	planetdata planet[MAXPLANETS];
	spotdata spot[MAXSPOTS];

	starttime=timecheck();

	CALCBRIGHTNESSFACTOR=1;	//unless changed in initializestarplanet
	FIXTHETAS=0;
	FIXSEEDEDONLYPHI=0;
	
	star=&thestar;
	
	if(argc<=1)
		sprintf(filename,DEFAULTFILENAME);
	else
		sprintf(filename,"%s",argv[1]);

	sprintf(rootname,"%s",filename);
	for(j=0;j<64&&rootname[j]!=0;j++);
	if(rootname[j-3]=='.'&&rootname[j-2]=='i'&&rootname[j-1]=='n')
		rootname[j-3]=0;

#	if DEBUGMCMC
		char dbfn[64];
		sprintf(dbfn,"%s_dbmcmc.txt",rootname);
		dbmcmc=fopen(dbfn,"w");
#	endif
#	if XYZDETAILS
		xyzdetail=fopen("xyzdetail.txt","w");
#	endif
		
	readparam=(double *)malloc(MAXSPOTS*3*sizeof(double));
	if(readparam==NULL) 
	{
		printf("malloc error\n");
		return 0;
	}
#	if !QUIET
		printf("initializing with parameters from %s\n",filename);
#	endif
	j=initializestarplanet(star,planet,filename,&lcstarttime,&lcfinishtime,&lcmaxlight,&ascale,&mcmcnpop,&mcmcmaxstepsortime,&partitionpop,&partitionsteps,readparam,&randomseed,seedfilename);
	srand(randomseed);
#	if PRINTEACHSPOT
		spotreport=(double *)malloc(numspots*sizeof(double));
		ldspotreport=(double *)malloc(numspots*sizeof(double));
		if(spotreport==NULL||ldspotreport==NULL)
			j=(-28);
#	endif
	if(j<0)
	{
		printf("initializestarplanet error %i\n",j);
		return 0;
	}

#	if !QUIET
		printf("\npreinitializing lc data from %s (t: %lf - %lf)\n",filename,lcstarttime,lcfinishtime);
#	endif
	preinitializelcdata(filename,lcstarttime,lcfinishtime,lcmaxlight,&lcn,&lclightnorm);
#	if !QUIET
		printf("lcn= %i\nlclightnorm= %lf\n",lcn,lclightnorm);
#	endif
	lctime=(double *)malloc(lcn*sizeof(double));
	lclight=(double *)malloc(lcn*sizeof(double));
	lcuncertainty=(double *)malloc(lcn*sizeof(double));
	if(lctime==NULL||lclight==NULL||lcuncertainty==NULL)
	{
		printf("malloc error\n");
		return 0;
	}
#	if !QUIET
		printf("initializing lc data\n");
#	endif
	i=initializelcdata(filename,lcstarttime,lcfinishtime,lcn,lclightnorm,lctime,lclight,lcuncertainty,planet);
	memused[0]=MAXSPOTS*(sizeof(tellipse)+sizeof(tcresent)+sizeof(tellipsehole)+sizeof(tcresenthole))+MAXPLANETS*sizeof(tcircle);
	memused[0]+=TWICEMAXSPOTSMAXPLANETS*(sizeof(tellipsesection)+sizeof(tcirclesection)+sizeof(tcresentsection)+sizeof(tellipsepart)+sizeof(tcresentpart)+sizeof(tpcirclearc)+sizeof(tlinesegment));
	memused[0]+=FOURMAXSPOTSMAXPLANETS*(sizeof(tellipsearc)+sizeof(tscirclearc));
	memused[1]=MAXSPOTS*3*sizeof(double)+3*lcn*sizeof(double);
#	if !QUIET
		printf("memory used for initialization: (%li+%li) = %li\n",memused[0],memused[1],memused[0]+memused[1]);
#	endif
	if(i<0)
	{
		printf("initialization error %i\n\n",i);
		return 0;
	}
#	if !QUIET
		else
			printf("initialization complete\n\n");
#	endif

#	if !QUIET
		for(i=0;i<2;i++)
			printf("%i: %lf %lf %lf\n",i,lctime[i],lclight[i],lcuncertainty[i]);
		printf("...\n");
		for(i=lcn-2;i<lcn;i++)
			printf("%i: %lf %lf %lf\n",i,lctime[i],lclight[i],lcuncertainty[i]);
#	endif

	sprintf(filename,"%s_errstsp.txt",rootname);
	outerr=fopen(filename,"w");
	if(outerr==NULL)
		printf("error opening outerr file: %s\n",filename);
	else
		fprintf(outerr,"--- error file for stsp ---\n");

#	if PRECALCPLANET
#		if !QUIET
			printf("initializing precalcplanet values\n");
#		endif
		for(i=0;i<MAXPLANETS;i++)
		{
			indcircle[i]=(tcircle *)malloc(lcn*sizeof(tcircle));
			if(indcircle[i]==NULL)
				printf("malloc error initializing indexed circle %i\n",i);
		}
		totalnumcircle=(int *)malloc(lcn*sizeof(int));
		if(totalnumcircle==NULL)
			printf("malloc error initializing totalnumcircle\n");
		initializeprecalcplanet(star,planet,lcn,lctime);
#		if !QUIET
			printf("  memory used: %li\n done initializing precalcplanet\n",lcn*sizeof(tcircle)*MAXPLANETS+lcn*sizeof(int));
#		endif
#	endif

	if(j==0||j==3||j==5||j==7||j==8||j==10)
	{
		PRINTVIS=0;
		if(j==0||j==8)
		{
#			if !QUIET
				printf("starting mcmc\n");
#			endif
			mcmc(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,ascale,mcmcnpop,mcmcmaxstepsortime,rootname,0,readparam,seedfilename);
		}
		else if(j==3||j==7)
		{
#			if !QUIET
				printf("starting seeded mcmc\n");
#			endif
			mcmc(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,ascale,mcmcnpop,mcmcmaxstepsortime,rootname,1,readparam,seedfilename);
		}
		else if(j==5)
		{
#			if !QUIET
				printf("starting totally seeded mcmc\n");
#			endif
			mcmc(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,ascale,mcmcnpop,mcmcmaxstepsortime,rootname,2,readparam,seedfilename);
		}
		else if(j==10)
		{
#			if !QUIET
				printf("starting partially seeded mcmc\n");
#			endif
			mcmc(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,ascale,mcmcnpop,mcmcmaxstepsortime,rootname,3,readparam,seedfilename);
		}
	}
	else if(j==1||j==6)
	{
#		if ANYPRINTVIS
		PRINTVIS=WHICHPRINTVIS;
		sprintf(filename,"%s_vis.txt",rootname);
		outv=fopen(filename,"w");
#		endif
		pvasc=1; pvisc=10;
		sprintf(seedfilename,"%s_lcout.txt",rootname);
		setspots(readparam,spot,star,planet);
		printf("generating light curve\n");
		lcgen(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,seedfilename);
	}
	else if(j==2)
	{
		PRINTVIS=0;
		printf("starting metropolis-hastings\n");
		metrohast(star,planet,spot,readparam,lcn,lctime,lclight,lcuncertainty,lclightnorm,mcmcmaxstepsortime,rootname);
	}
	else if(j==4)
	{
		printf("starting debug trials\n");
		debuglc(star,planet,spot,lcn,lctime,lclight,lcuncertainty);
	}
	else if(j==9)
	{
		printf("starting plot data for one spot variation\n");
		plotdata(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,randomseed,readparam,rootname);
	}
	else if(j==11)
	{
		initialguess(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,seedfilename);
	}
	else if(j==12||j==13)
	{
		sprintf(seedfilename,"%s_acout.txt",rootname);
		setspots(readparam,spot,star,planet);
		printf("finding spot area not covered\n");
		i=areacovered(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,seedfilename);
		if(i)
			printf("error %i\n",i);
	}
	else if(j==100)
	{
		gentimetest(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,0,readparam);
	}
	else if(j==101)
	{
		gentimetest(star,planet,spot,lcn,lctime,lclight,lcuncertainty,lclightnorm,1,readparam);
	}

	fclose(outerr);
	free((void *)readparam);
	free((void *)lctime);
	free((void *)lclight);
	free((void *)lcuncertainty);
#	if PRINTEACHSPOT
		free((void *)spotreport);
		free((void *)ldspotreport);
#	endif
#	if DEBUGMCMC
		fclose(dbmcmc);
#	endif
#	if XYZDETAILS
		fclose(xyzdetail);
#	endif

	return 0;
}


