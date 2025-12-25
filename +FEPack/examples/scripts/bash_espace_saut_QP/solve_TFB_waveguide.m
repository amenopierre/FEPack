function solve_TFB_waveguide(idBEGIN, idEND, numNodes, cheminDonnees)

  % Load inputs
  load([cheminDonnees, '/inputs_', int2str(numNodes), '.mat']);

  % Restrict idEND to the number of Floquet points
  idEND = min(idEND, numFloquetPoints);
  
  % Solve the waveguide problems
  for idFB = idBEGIN:idEND
    % tic;
    fprintf('%d sur %d\n', idFB, idEND-idBEGIN+1);
    
    FloquetVar = FloquetPoints(idFB);

    % The FE matrix is a linear combination of the elementary pieces
    fprintf('Calcul matrices EF\n');
    AApos = mat_gradu_gradv_pos + 1i * FloquetVar * mat_vec1u_gradv_pos - 1i * FloquetVar * mat_gradu_vec1v_pos + FloquetVar * FloquetVar * mat_vec1u_vec1v_pos - (opts.omega * opts.omega) * mat_u_v_pos;
    AAneg = mat_gradu_gradv_neg + 1i * FloquetVar * mat_vec1u_gradv_neg - 1i * FloquetVar * mat_gradu_vec1v_neg + FloquetVar * FloquetVar * mat_vec1u_vec1v_neg - (opts.omega * opts.omega) * mat_u_v_neg;

    % The Floquet-Bloch transform of the boundary data
    fprintf('Calcul TFB donnee de saut\n');
    jumpData_FB = @(x) FEPack.tools.BlochTransform(x, FloquetVar, G3D, infiniteDirection);

    % Compute the Floquet-Bloch transform of the solution
    fprintf('Calcul solution probleme de transmission\n');
    TFBU = PeriodicGuideJumpBVP(semiInfiniteDirection,...
                                AApos, mesh3Dpos, BCstruct_pos, numCellsSemiInfinite_pos,...
                                AAneg, mesh3Dneg, BCstruct_neg, numCellsSemiInfinite_neg,...
                                jumpData_FB, opts);
    
    % Save the Floquet-Bloch transform of the solution
    save([cheminDonnees, '/TFBU_', int2str(idFB), '.mat'], '-struct', 'TFBU', '-v7.3');
    % toc;
  end
  
end