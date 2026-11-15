// Переходный канал (или улитка) KV: adjointOptimisationFoam openfoam.com v2512 - шаблон optimisationDict.
// Объекты и ключи - как в учебниках v2512 sbend (PtLosses, volumetricBSplines, BFGS) и fork-uneven (uniformityPatch).
optimisationManager     steadyOptimisation;
primalSolvers
{
    p1 { active true; type incompressible; solver simple;
         solutionControls { nIters 3000; residualControl { "p.*" 1e-7; "U.*" 1e-7; } } }
}
adjointManagers
{
    am1
    {
        primalSolver p1;
        adjointSolvers
        {
            as1
            {
                active true; type incompressible; solver adjointSimple;
                objectives
                {
                    type incompressible;
                    objectiveNames
                    {
                        losses     { weight 1;   type PtLosses;        patches (inlet outlet); }
                        uniformity { weight <w>; type uniformityPatch; patches (outlet); }
                    }
                }
                ATCModel { ATCModel standard; }
                solutionControls { nIters 3000; residualControl { "pa.*" 1e-7; "Ua.*" 1e-7; } }
            }
        }
    }
}
optimisation
{
    designVariables { type shape; shapeType volumetricBSplines; sensitivityType shapeFI; patches (walls); maxInitChange 2e-3; }
    updateMethod { method BFGS; BFGS { etaHessian 0.8; nSteepestDescent 1; scaleFirstHessian true; } }
}
// constant/dynamicMeshDict: solver volumetricBSplinesMotionSolver; box: type cylindrical; nCPsU <nu>; nCPsV <nv>; nCPsW <nw>
