DefineConstant[
    mu0 = 4*Pi*1e-7,
    DSVsel = {DSV, Min 0.001, Max 1, Step 0.001,
        Name Sprintf("Input/0DSV")},
    Flag_Analysis = {0,
        Choices { 0 = "H-formulation", 1 = "A-formulation"},
        Name "Input/1Formulation" },
    NumMagnets_ = {NumMagnets,
        Name Sprintf("Input/5Number of magnets"),
        ReadOnly 1}    
];

Group{
  Vol_Magnets_Mag = Region[{}];
  For i In {1:NumMagnets}
    Magnet~{i} = Region[i]; // volume of magnet i
    SkinMagnet~{i} = Region[(SurfaceRegionOffset+i)]; // boundary of magnet i
    Vol_Magnets_Mag += Region[Magnet~{i}]; // all the magnet volumes
  EndFor
  Air = Region[(NumMagnets + 1)];
  Vol_DSV = Region[(NumMagnets + 3)];
  Vol_Mag = Region[{Air, Vol_Magnets_Mag}];
  Dirichlet_phi_0 = Region[(NumMagnets + 2)]; // boundary of air box
  Dirichlet_a_0 = Region[(NumMagnets + 2)]; // boundary of air box
}

Function{
  mu[Air] = mu0;
  nu[Air] = 1.0/mu0;
  mu[Vol_DSV] = mu0;
  nu[Vol_DSV] = 1.0/mu0;
  For i In {1:NumMagnets}
    // coercive field of magnets

    hc[Magnet~{i}] = Rotate[Vector[BR~{i}/mu0, 0, 0], 0, 0, angle~{i}];
    br[Magnet~{i}] = Rotate[Vector[BR~{i}, 0, 0], 0, 0, angle~{i}];
    mu[Magnet~{i}] = mur~{i} * mu0;
    nu[Magnet~{i}] = 1.0/(mur~{i} * mu0);
  EndFor
}

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
	  { GeoElement Triangle ; NumberOfPoints 4 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
          { GeoElement Tetrahedron  ; NumberOfPoints 4 ; }
	}
      }
    }
  }
}

Constraint {
  { Name phi ;
    Case {
      { Region Dirichlet_phi_0 ; Value 0. ; }
    }
  }
  { Name a ;
    Case {
      { Region Dirichlet_a_0 ; Value 0. ; }
    }
  }
  { Name GaugeCondition_a ; Type Assign ;
    Case {
      { Region Vol_Mag ; SubRegion Dirichlet_a_0 ; Value 0. ; }
    }
  }
}

FunctionSpace {
  // scalar magnetic potential
  { Name Hgrad_phi ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef phin ; Function BF_Node ;
        Support Vol_Mag ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef phin ; EntityType NodesOf ; NameOfConstraint phi ; }
    }
  }
  // vector magnetic potential
  { Name Hcurl_a; Type Form1;
    BasisFunction {
      { Name se;  NameOfCoef ae;  Function BF_Edge; Support Vol_Mag ;
        Entity EdgesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae;  EntityType EdgesOf ; NameOfConstraint a; }
      { NameOfCoef ae;  EntityType EdgesOfTreeIn ; EntitySubType StartingOn ;
        NameOfConstraint GaugeCondition_a ; }
    }
  }
}

Formulation {
  { Name MagSta_phi ; Type FemEquation ;
    Quantity {
      { Name phi ; Type Local ; NameOfSpace Hgrad_phi ; }
    }
    Equation {
      Integral { [ - mu[] * Dof{d phi} , {d phi} ] ;
        In Vol_Mag ; Jacobian JVol ; Integration I1 ; }
      Integral { [ - mu[] * hc[] , {d phi} ] ;
        In Vol_Magnets_Mag ; Jacobian JVol ; Integration I1 ; }
    }
  }
  { Name MagSta_a; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a ; }
    }
    Equation {
      Integral { [ nu[] * Dof{d a} , {d a} ] ;
        In Vol_Mag ; Jacobian JVol ; Integration I1 ; }
      Integral { [ nu[] * br[] , {d a} ] ;
        In Vol_Magnets_Mag ; Jacobian JVol ; Integration I1 ; }
    }
  }
}

Resolution {
  { Name MagSta_phi ;
    System {
      { Name A ; NameOfFormulation MagSta_phi ; }
    }
    Operation {
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
      PostOperation[MagSta_phi] ;
    }
  }
  { Name MagSta_a ;
    System {
      { Name A ; NameOfFormulation MagSta_a ; }
    }
    Operation {
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
      PostOperation[MagSta_a] ;
    }
  }
  { Name Analysis ;
    System {
      If(Flag_Analysis == 0)
        { Name A ; NameOfFormulation MagSta_phi ; }
      EndIf
      If(Flag_Analysis == 1)
        { Name A ; NameOfFormulation MagSta_a ; }
      EndIf
    }
    Operation {
      Generate[A] ; Solve[A] ; SaveSolution[A] ;
      If(Flag_Analysis == 0)
        PostOperation[MagSta_phi] ;
      EndIf
      If(Flag_Analysis == 1)
        PostOperation[MagSta_a] ;
      EndIf
    }
  }
}

PostProcessing {
  { Name MagSta_phi ; NameOfFormulation MagSta_phi ;
    PostQuantity {
      { Name b ; Value {
          Term { [ - mu[] * {d phi} ] ; In Vol_Mag ; Jacobian JVol ; }
          Term { [ - mu[] * hc[] ]    ; In Vol_Magnets_Mag ; Jacobian JVol ; }
        }
      }
      { Name b_abs ; Value {
          Term { [Norm[ - mu[] * {d phi}]]  ; In Vol_Mag ; Jacobian JVol ; }
        }
      }
      { Name b_x ; Value {
          Term { [Fabs[CompX[ - mu[] * {d phi}]]]  ; In Vol_Mag ; Jacobian JVol ; }
        }
      }
      { Name h ; Value {
          Term { [ - {d phi} ] ; In Vol_Mag ; Jacobian JVol ; }
        }
      }
      { Name hc ; Value {
          Term { [ hc[] ] ; In Vol_Magnets_Mag ; Jacobian JVol ; }
        }
      }
      { Name phi ; Value {
          Term { [ {phi} ] ; In Vol_Mag ; Jacobian JVol ; }
        }
      }
    }
  }
  { Name MagSta_a ; NameOfFormulation MagSta_a ;
    PostQuantity {
      { Name b ; Value {
          Term { [ {d a} ]; In Vol_Mag ; Jacobian JVol; }
        }
      }
      { Name a ; Value {
          Term { [ {a} ]; In Vol_Mag ; Jacobian JVol; }
        }
      }
      { Name br ; Value {
          Term { [ br[] ]; In Vol_Magnets_Mag ; Jacobian JVol; }
        }
      }
    }
  }
}

PostOperation {
  { Name MagSta_phi ; NameOfPostProcessing MagSta_phi;
    Operation {
      Print[ b, OnElementsOf Vol_Mag, File StrCat[outputFilename, "_b.pos"] ] ;
      Print[ b, OnPlane{ {-DSVsel/2, -DSVsel/2, 0} {DSVsel/2, -DSVsel/2, 0} {-DSVsel/2, DSVsel/2, 0} } {200, 200},
        File "b_cut1.pos" ];
      Print[ b_abs, OnGrid {$B*Cos[$A], $B*Sin[$A], 0} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_abs_frontal.pos"]];
      Print[ b_x, OnGrid {$B*Cos[$A], $B*Sin[$A], 0} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_x_frontal.pos"]];
      Print[ b_abs, OnGrid {$B*Cos[$A], 0, $B*Sin[$A]} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_abs_transverse.pos"]];
      Print[ b_x, OnGrid {$B*Cos[$A], 0, $B*Sin[$A]} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_x_transverse.pos"]];
      Print[ b_abs, OnGrid {0, $B*Cos[$A], $B*Sin[$A]} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_abs_longitudinal.pos"]];
      Print[ b_x, OnGrid {0, $B*Cos[$A], $B*Sin[$A]} { 0:2*Pi:2*Pi/1000, 0:DSVsel/2:DSVsel/2/100, 0 },
        File StrCat[outputFilename, "_B_x_longitudinal.pos"]];
      //Print[ h, OnElementsOf Vol_Mag, File "h.pos" ] ;
      //Print[ hc, OnElementsOf Vol_Mag, File "hc.pos" ] ;
    }
  }
  { Name MagSta_a ; NameOfPostProcessing MagSta_a ;
    Operation {
      Print[ b,  OnElementsOf Vol_Mag,  File "b.pos" ];
      Print[ b, OnPlane{ {-0.1,-0.1,0} {0.1,-0.1,0} {-0.1,0.1,0} } {50, 50},
        File "b_cut1.pos" ];
      //Print[ br,  OnElementsOf Vol_Magnets_Mag,  File "br.pos" ];
      //Print[ a,  OnElementsOf Vol_Mag,  File "a.pos" ];
    }
  }
}

DefineConstant[
  // preset all getdp options and make them invisible
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -v 5 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 0}
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];
