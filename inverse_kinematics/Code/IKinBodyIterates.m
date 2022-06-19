%this function: IKinBodyIterates(Blist, M, T, thetalist0, eomg, ev) takes
% the following inputs:
% 1) Blist: the matrix containing the screw axises expressed in the body
% frame as its columns.
% 2) M: the zero-position home configuration of end-effector expressed in
% the space frame.
% 3) T: the desaired end-effector config. expressed in the space frame.
% 4) thetalist0: the initial guess column vector of joint variables.
% 5) eomg and ev: the tolerance in the config. of the resulted end
% effector frame. it's represented as the error of the body twist of the end
% effector. The magnitude of both the omega and v vectors is zero at the
% exact config. desaired.
% the output of each iteration is: the joint variables vector, the config,
% the body twist and the errors in omega and v.
% the final outputs are the solution, the joint vector, and a matrix
% containing all previous iterations joint vectors as its rows and saved as
% a file named iterates.csv at your current MATLAB workplace directory.
function [thetalist, success] = IKinBodyIterates(Blist, M, T, thetalist0, eomg, ev)
  thetalist = thetalist0; % this command initializes the variable thetalist given the initial guess..
  n=size(thetalist0,1);   % n is a variable returns the first dimension of the vector thetalist0 which is the joints number.
  maxiterations = 20;     
  thetaMat=zeros(maxiterations+1,n); % initializing a thetaMat matrix of dimension maxiteration+1 for rows, and n for cols. 
  thetaMat(1,:)=transpose(thetalist0); % the first column of thetaMat contains the initial guess. 
  i = 0;
  Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T)); % calculating the error twist Vb for the initial guess.
  err = norm(Vb(1: 3)) > eomg || norm(Vb(4: 6)) > ev; % calc the error in omega and v magnitudes giving the initial guess.
  while err && i < maxiterations % iterate the following while the error is larger than the specified input value.
      thetalist = thetalist + pinv(JacobianBody(Blist, thetalist)) * Vb;
      i = i + 1;
      Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T));
      magWb = norm(Vb(1: 3));  % calculate the new error in omega.
      magVb = norm(Vb(4: 6));  % calculate the new error in vb.
      % calc the current matrix of the end effector config respective to space frame.
      Ti=FKinBody(M, Blist, thetalist); 
      % print out each iteration number and the corresponding joint variables.
      fprintf('Iteration:\n%1d\nJoint vector:\n%4f %4f %4f %4f %4f %4f\n',i,thetalist); 
      % display related configuration.
      fprintf('SE(3) end effector configuration(with respect to space frame): \n \n');
      disp(Ti); 
      % display the related error twist vb,and the related errors.
      Text2='The error twist Vb:\n%4f %4f %4f %4f %4f %4f\nAngular error magnitude:\n%4f\nLinear error magnitude:\n%4f\n';  
      fprintf(Text2,Vb,magWb,magVb);
      disp('---------------------------------------------------------------');
      err= magWb > eomg || magVb > ev; % check whether the input tolerance in the error is reached.
      thetaMat(i+1,:)=transpose(thetalist); % fill the thetaMat rows by the current joint vector thetalist.
  end
  thetaMat(i+2:end,:) = []; % remove the remaining unused rows in thetaMat starting from (i+2th row.
  % print out the joint iteration matrix and its size.
  fprintf('The joint iterations matrix:\n\n'); disp(thetaMat);
  fprintf('Joint matrix size is %4d Rows(representing iterations) X %4d Cols(representing joint variables). \n',size(thetaMat));
  % finally, save the iteration matrix (named iterates.csv) in the current MATLAB workspace
  csvwrite('iterates.csv',thetaMat);
  success = ~ err;
end
