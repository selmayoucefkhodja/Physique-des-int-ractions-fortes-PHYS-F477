#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TRandom3.h>
#include <time.h>
#include <cmath>

// On definit des constantes valables
// dans tout le code.
#define pi 3.141592
#define Ebeam 10.
#define alphaEM 0.0072973525
#define alphaS 0.1184
#define Nc 3
#define Cf 4./3.
#define s 4*Ebeam*Ebeam
#define Qq 2./3.

using namespace std;

// Codez ici une formule permettant de faire le passage
// du referentiel alpha, beta, gamma au referentiel du
// laboratoire.
// Cette fonction doit prendre les 3 angles en argument
// plus une reference au quadrivecteur a modifier.
void myRotate(double alpha, double beta, double gamma, TLorentzVector &particle)
{
  // Pour acceder au trivecteur du quadrivecteur
  // on utilise la methode Vect()
  TVector3 myVect = particle.Vect();

  TVector3 X (1,0,0);                        //Définition des axes de base dans le référentiel d'Euler.
  TVector3 Z (0,0,1);

  myVect.Rotate(-gamma,Z);                  //Le signe négatif est lié aux définitons des angles alpha beta gamma pour le passage du référentiel du labo au régérentiel d'Euler (là, on fait le chemin inverse).
  X.Rotate(-gamma,Z);                       //On fait tourner les axes du référentiel d'Euler puisque les rotations sont définies par rapport au référentiel obtenu après chaque étape.
  myVect.Rotate(-beta,X);
  Z.Rotate(-beta, X);
  myVect.Rotate(-alpha,Z);

  // La methode Rotate(alpha, V) effectue
  // une rotation d un angle alpha autour de
  // l axe defini par le vecteur V

  particle.SetVect(myVect);
}

// Libre a vous de tester la fonction
// de rotation que vous venez de coder

void testMyRotate()
{
  TLorentzVector ele(1.,0.,0.,0.);  //(px, py, pz, E) dans le ref. alpha, beta, gamma
  myRotate(pi/4, pi/4, 0, ele);
  cout << ele.Px() << " " << ele.Py() << " " << ele.Pz() << endl;
  TLorentzVector pos(0.,1.,0.,0.);
  myRotate(pi/4, 0, 0, pos);
  cout << pos.Px() << " " << pos.Py() << " " << pos.Pz() << endl;
  TLorentzVector posbis(1.,1.,1.,0.);
  myRotate(pi/10, pi/6, pi/4, posbis);
  cout << posbis.Px() << " " << posbis.Py() << " " << posbis.Pz() << endl;
}

// Codez ici l expression de la section efficace
// differentielle en termes des quadrivecteurs.
  double dsigma(TLorentzVector ele, TLorentzVector pos, TLorentzVector quark, TLorentzVector antiq, TLorentzVector gluon)
{
  double result = (Cf*alphaEM*alphaEM*alphaS*Nc*Qq*Qq/(2*pi*pi))*( (pos*quark)*(pos*quark) +(pos*antiq)*(pos*antiq) + (ele*antiq)*(ele*antiq) )/(s*(quark*gluon)*(antiq*gluon));   //Juste la formule du TP4
  //result = factor * bracket / deno;
  return result;
}

// Codez ici l expression de la section efficace
// differentielle en termes des variables x1, x2, beta, gamma
  double dsigma(double x1, double x2, double beta, double gamma)
{
  double Theta1, Theta2, x3;
  x3 = 2. - x1 - x2;
  Theta1 = acos(-(x1*x1 + x3*x3 - x2*x2)/(2.*x1*x3));
  Theta2 = acos(-(x2*x2 + x3*x3 - x1*x1)/(2.*x2*x3));

  double factor = Cf*alphaEM*alphaEM*alphaS*Nc*Qq*Qq/(4*pi*pi);
  double bracket = (x1*x1*(1+sin(beta)*sin(beta)*sin(gamma + Theta1)*sin(gamma + Theta1))+x2*x2*(1+sin(beta)*sin(beta)*sin(gamma - Theta2)*sin(gamma - Theta2)));
  double denom = (s*(1-x1)*(1-x2));

  double result =  factor * bracket / denom;
  //result = factor*bracket/deno;
   return result;
}
void testDsigma()
{
   cout << dsigma(0.8, 0.8, pi/4, pi/4) << endl;
}

// Codez ici la fonction principale qui genere les
// evenements et qui applique la methode de rejection.
void rejMethode()
{
  // Determinez le sigmaMax en faisant une premiere
  // boucle et en sauvant la valeur maximale de la
  // section efficace
  double sigmaMax =  2.83E-9;                                                           //Initialisation du maximum de la séction efficace.

  // Declarez toutes les variables necessaires pour
  // votre programme
  int npoints = 50000000, nAccepted = 0;
  double alpha, cosBeta, beta, gamma, x1, x2, x3, cosTheta1, cosTheta2, fmax, y1,y2,y_min,y_max;
  double x1_min = 0.00005;
  double number1jetevents = 0;                                                         // Pour vérifier que l'on n'a pas d'événements à un seul jet.

  // Declarez les 5 quadri vecteurs ainsi que les
  // variables a sauver
  TLorentzVector ele, pos, quark, antiq, gluon;
  double x1_Rec, x2_Rec, x3_Rec;



  // Declarez une serie d histogrammes vous permettant
  // de controler vos resultats
  TH1F *h_cosThetaQuark = new TH1F("h_cosThetaQuark","h_cosThetaQuark;cos#theta_{quark}",100, -1., 1.);           //Histogrammes des variables \theta_{1,2} après application de la méthode de réjection.
  TH1F *h_cosThetaAntiq = new TH1F("h_cosThetaAntiq","h_cosThetaAntiq;cos#theta_{antiquark}",100, -1., 1.);
  TH1F *h_alpha = new TH1F("h_alpha","h_alpha;#alpha",100, 0., 6.29);                  //Histogrammes des variables angulaires d'Euler avant application de la méthode de réjection.
  TH1F *h_gamma = new TH1F("h_gamma","h_gamma;#gamma",100, 0., 6.29);
  TH1F *h_cosBeta = new TH1F("h_cosBeta","h_cosBeta;cos#beta",100, -1., 1.);
  TH1F *h_beta = new TH1F("h_beta","h_beta;#beta",100, 0., 3.142);

  TH1F *h_alphar = new TH1F("h_alphar","h_alpha;#alpha",100, 0., 6.29);                 //Histogrammes des variables angulaires d'Euler après application de la méthode de réjection.
  TH1F *h_gammar = new TH1F("h_gammar","h_gamma;#gamma",100, 0., 6.29);
  TH1F *h_cosBetar = new TH1F("h_cosBeta","h_cosBetar;cos#beta",100, -1., 1.);
  TH1F *h_betar = new TH1F("h_betar","h_beta;#beta",100, 0., 3.142);

  TH3F *h_xxx = new TH3F("h_xxx", "h_xxx;x_1;x_2;x_3", 100, 0, 1, 100, 0, 1, 100, 0, 1);//Corrélation entre les x_i ainsi que leurs distributions respectives.
  TH1F *h_x1 = new TH1F("h_x1","h_x1;x_1",100,0.5,1);
  TH1F *h_x2 = new TH1F("h_x2","h_x2;x_2",100,0.4,1);
  TH1F *h_x3 = new TH1F("h_x3","h_x3;x_3",100,0,1);

  TH1F *h_x1b = new TH1F("h_x1b","h_x1;x_1",100,0.5,1);                                 //Hisogrammes des x_1 générées avant l'application de la méthode de réjection. Utile pour comprendre l'évolution de l'efficacité du code.

  TH1F *h_x1_r = new TH1F("h_x1r","h_x1r;x_1r",100,0.5,1);                              //Histogrammes de x_i avant reconstruction. Utile pour comparer.
  TH1F *h_x2_r = new TH1F("h_x2r","h_x2r;x_2r",100,0.4,1);
  TH1F *h_x3_r = new TH1F("h_x3r","h_x3r;x_3r",100,0,1);
  TH2F *h_x1x2 = new TH2F("h_x1x2","h_x1x2;x_{1};x_{2}",100,0,1,100,0,1);
  TH1F *h_cosThetaEK = new TH1F("h_cosThetaEK","h_cosThetaEK;cos#theta_{EK}",100,0,1);  //Histogramme de cos\theta_{EK} où \theta_{EK} correspond à l'angle entre les jets les plus énergetiques.
  TH2F *h_sigmax1 = new TH2F("h_sigmax1","h_sigmax1;#sigma;x_1",1000, 0, 0.000000004, 100, 0, 1);    //Histogramme avec x1 après reconstruction.
  TH2F *h_sigmax1r = new TH2F("h_sigmax1r","h_sigmax1;#sigma;x_1",1000, 0, 0.000000004, 100, 0, 1);  //Histogramme avec x1 avant reconstruction. Utile pour comparer.
  TH1F *h_2jets = new TH1F("h_2jets","h_2jets;y_{cut}",100, 0., 1.);                    //Comme on s'intérésse aux événements e+ e- -> q qbar g, on a maximum 3 jets et minimum 2 jets.
  TH1F *h_3jets = new TH1F("h_3jets","h_3jets;y_{cut}",100, 0., 1.);
  // Initialisez le generateur de nombres aleatoires
  TRandom3 generateur;
  generateur.SetSeed(123456);

  // Ici commence la boucle sur les evenements a generer
  for (int i=0; i<npoints; i++)
  {
   // Generez les differents nombres aléatoires necessaires
   // sur les intervalles adequats

      //x1 = (1-2*x1_min)*generateur.Rndm() + x1_min ;           //Ancienne façon pas efficace de générer x1 et x2. Or la section efficace contribue bien plus pour des valeurs proches de 1 car on a une divergence. Donc il vaut
                                                                 //mieux génerer plus d'événements dans cette région là afin d'augmenter l'efficacité du code.
      //x2 = (1-2*x1_min)*generateur.Rndm() + x1_min ;

      y_min = -log(1-x1_min);                                    //On utilise un changement de variable afin d'améliorer l'efficacité du code, il fait en sorte que la distribution de la section efficace est plus applatie.
      y_max = - log(x1_min);                                     //Ce changement de variable impose également de changer les intervalles parcourus.
      y1 = (y_max - y_min)*generateur.Rndm() + y_min;
      y2 = (y_max - y_min)*generateur.Rndm() + y_min;
      x1 = 1.- exp(-y1);                                         //On inverse la relation du changement de variable afin d'avoir x en fonction de y.
      x2 = 1.- exp(-y2);
      x3 = 2.- x1 - x2;
      if(x3 >= 1.-x1_min) continue;                               //Car si x1 et x2 prennent leur valeur minimale x1_min, x3 = 2 -2*x2_min alors qu'il ne doit pas dépasser 1-x1_min. Par contre, si x1 et x2 prennent leur valeur
                                                                  //maximale 1-x1_min, on trouve un x3 prend la valeur minimale  2 - (1- x1_min) - (1 - x1_min) = 2*x1_min qui est bien supérieure à x_min. Donc pas besoin d'imposer de conditions supplémentaires.
      alpha = (2*pi)*generateur.Rndm() ;
      cosBeta = 2*generateur.Rndm() - 1 ;
      beta = acos(cosBeta);
      gamma = 2*pi*generateur.Rndm() ;                            //Formules et définitions.
      cosTheta1 = -(x1*x1 + x3*x3 - x2*x2)/(2*x1*x3);
      cosTheta2 = -(x2*x2 + x3*x3 - x1*x1)/(2*x2*x3);
      fmax = sigmaMax*generateur.Rndm();                          //utile pour la méthode de rejection.



   // Construisez toute la cinematique dans le referentiel alpha, beta, gamma
   TLorentzVector ele(Ebeam*sqrt(1-cosBeta*cosBeta)*sin(gamma),Ebeam*sqrt(1-cosBeta*cosBeta)*cos(gamma),Ebeam*cosBeta,Ebeam);     //résultats du TP 4. Cela dit, seuls les quadrivecteurs des quark/antiquark/gluon interviendront.
   TLorentzVector pos(-Ebeam*sqrt(1-cosBeta*cosBeta)*sin(gamma),-Ebeam*sqrt(1-cosBeta*cosBeta)*cos(gamma),-Ebeam*cosBeta,Ebeam);
   TLorentzVector quark( Ebeam*x1*cosTheta1,Ebeam*x1*sqrt(1-cosTheta1*cosTheta1),0,Ebeam*x1*1);
   TLorentzVector antiq(Ebeam*x2*cosTheta2,-Ebeam*x2*sqrt(1-cosTheta2*cosTheta2),0,Ebeam*x2*1);
   TLorentzVector gluon(Ebeam*(2-x1-x2),0,0,Ebeam*(2-x1-x2));


    // Passage au referentiel du labo
    myRotate(alpha, beta, gamma, quark);
    myRotate(alpha, beta, gamma, antiq);
    myRotate(alpha, beta, gamma, gluon);

    h_alpha->Fill(alpha); h_gamma->Fill(gamma);	h_beta->Fill(beta); h_cosBeta->Fill(cosBeta); //Avant l'application de la méthode de réjection. Utile pour comparer
    h_x1b -> Fill(x1);

    // Inserez la condition de la methode de rejection

    double Sigma = dsigma(x1,x2,beta,gamma)*(1-x1)*(1-x2);   //On multiplie par le jacobien pour calculer la section efficace avec le nouveau changement de variable.
    if (Sigma > sigmaMax) {                                  //Pour mettre à jour et afficher  la valeur maximale de la section efficace differentielle
      sigmaMax = Sigma;
      cout << sigmaMax  << endl;
    }
    if ( Sigma < fmax ) continue;                            //La condition de rejection nous dit qu'on doit rejeter l'événement si le fmax généré uniformément dépasse sigma évaluée en les autres variables générées uniformément.
    nAccepted++;

     x1_Rec = max( max(x1,x2),x3 );                          //Comme le détécteur ne sait pas quelle fraction est prise par quelle particule, on définit ces nouvelles variables rangées par ordre de croissance.
     x3_Rec =  min( min(x1,x2),x3 );
     x2_Rec = 2 - x3_Rec - x1_Rec;

     // Remplissez vos Histogrammes
  h_alphar->Fill(alpha); h_gammar->Fill(gamma);	h_betar->Fill(beta); h_cosBetar->Fill(cosBeta);   //Après l'application de la méthode de réjection.

  h_xxx -> Fill(x1_Rec,x2_Rec, x3_Rec);                       //Correlation entre les trois variables x1,x2,x3.
  h_x1 -> Fill(x1_Rec);                                       //On remplit les histogrammes des x_i après reconstruction.
  h_x2 -> Fill(x2_Rec);
  h_x3  -> Fill(x3_Rec);

  h_x1_r -> Fill(x1);                                         //On remplit les histogrammes des x_i avant reconstruction. Utile pour comparer.
  h_x2_r -> Fill(x2);
  h_x3_r  -> Fill(x3);
  h_x1x2  -> Fill(x1_Rec,x2_Rec);
  h_cosThetaEK -> Fill((x2_Rec-x3_Rec)/x1_Rec);
  h_sigmax1 -> Fill(Sigma, x1_Rec);                            //x1 reconstruit comme fonction de la section efficace.
  h_sigmax1r -> Fill(Sigma, x1);                               //x1 avant reconstruction comme fonction de la section efficace. Utile pour comparer.
  h_cosThetaQuark->Fill(cos(quark.Theta())); h_cosThetaAntiq->Fill(cos(antiq.Theta()));


   //Algorithme de Jade
  double ycut, jetnumber;
  TLorentzVector p1, p2,p3,p4;
  jetnumber = 3;                                              //Car on ne peut pas avoir plus de 3 jets dans ce contexte.
  ycut = generateur.Rndm();                                   //Car on veut des histogrammes de distribution des évènements selon le ycut choisi.
  p1 = quark + antiq;                                         //Définitions utile pour tester le critère d'appartenance à un même jet.
  p2 = quark + gluon;
  p3 = antiq + gluon;
  p4 = quark + antiq + gluon;
                                                            //Si une des conditions suivantes est respectée, on a forcement 2 jets pcq il y aura une pair de particules (soit q qbar ou q gluon ou qbar gluon) qui vont former une meta particule et donc le jetnumber baissera d'une unité.
                                                            //De plus, on sait qu'on ne peut pas descendre plus bas parce que on aura jamais p4*p4<ycut*s parce que p4*p4= s par conservation de l'energie et au mieux ycut=1 donc la condition s<s ne sera jamais satisfaite.
  if(p1*p1 < ycut*s or p2*p2< ycut*s or p3*p3 < ycut*s ){
  h_2jets -> Fill(ycut);
}
  else{                                                      //Logiquement, si toutes les paires de particules ont une énergie supérieur au threshold imposé, on ne "fusionnera" aucune paire et on aura bien 3 jets.
  h_3jets -> Fill(ycut);
}

  if (p4*p4 < ycut*s){                                       // Juste pour bien voir qu'il n'y a aucun événement à un seul jet !
  number1jetevents = number1jetevents + 1;
}


}
  //fin de la boucle de generation d evenements
  cout << "Nombre d'évnement à un seul Jet = " << number1jetevents << endl;            // Pour bien voir que l'on a aucun événement à un seul jet.

  // Hors de la boucle, dessinez vos histogrammes dans
  // des canevas. Faites des fit sur certaines distributions,
  // sauvez les histogrammes.

  TF1 *f = new TF1("f","[0]*(1.+[1]*x*x)", -0.9, 0.9);                                 //Fit avec une fonction bien symétrique.
  f->SetParameters(100, 1);
  TF1 *g = new TF1("g","[0*x]/(1.-[1]*x)", 1, 1);                                     //Fit avec une fonction en x/(1-x).
  g->SetParameters(1, 1);

  TCanvas *can = new TCanvas(); can->Divide(2,4);
  can->cd(1);	h_cosThetaQuark->Fit("f"); h_cosThetaQuark->DrawCopy();              //Histogrammes des variables angulaires \theta_{EK}, \alpha, \beta, \gamma, et \cos\beta.
  can->cd(2);	h_cosThetaAntiq->Fit("f"); h_cosThetaAntiq->DrawCopy();
  can->cd(3);	h_alpha->DrawCopy();
  can->cd(4);	h_beta->DrawCopy();
  can->cd(5);	h_gamma->DrawCopy();
  can->cd(6);	h_cosBeta->DrawCopy();
  can->cd(7);	h_2jets->DrawCopy();                                                //Histogramme de la distribution des événements à 2/3 jets en fonction du y_{cut}.
  can->cd(8);	h_3jets->DrawCopy();

  TCanvas *can2 = new TCanvas(); can2->Divide(2,4);                                 //Pour afficher les histogrammes faisant intervenir les fraction d'énergie reconstruites.
  can2->cd(1); h_xxx->DrawCopy();
  can2->cd(2); h_x1->DrawCopy();
  can2->cd(3); h_x2->DrawCopy();
  can2->cd(4); h_x3->DrawCopy();
  can2->cd(5); h_sigmax1->DrawCopy();
  can2->cd(6); h_x1x2->DrawCopy();
  can2->cd(7); h_cosThetaEK->DrawCopy();
  can2->cd(8); h_x1b->DrawCopy();

  TCanvas *can3 = new TCanvas(); can3->Divide(2,4);
  can3->cd(1); h_x1_r->DrawCopy()->Fit("g");                                         //Pour afficher les histogrammes faisant intervenir les fractions d'énergie non reconstruites.
  can3->cd(2); h_x2_r->DrawCopy();
  can3->cd(3); h_x3_r->DrawCopy();
  can3->cd(4);	h_alphar->DrawCopy();                                                //Pour afficher la distribution des angles \alpha, \beta, \gamma et \cos\beta après application de la méthode de réjection.
  can3->cd(5);	h_betar->DrawCopy();
  can3->cd(6);	h_gammar->DrawCopy();
  can3->cd(7);	h_cosBetar->DrawCopy();
  can3->cd(8); h_sigmax1r->DrawCopy();                                               //La fraction d'énergie x_1 non reconstruite comme fonction de la section efficace.


  cout << "Sigma Max = " << sigmaMax << endl;
  cout << "Eff = " << nAccepted*1./npoints << endl;

}
