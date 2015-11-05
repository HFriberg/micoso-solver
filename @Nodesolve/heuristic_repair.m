function heuristic_repair(self, prob)
  %
  % Currently only repairing in primal variables
  %
  
  % Assuming no more than one slack per row!
  slackcol = find(prob.c(1:prob.K.l)==0 & sum(prob.A(:, 1:prob.K.l) ~= 0)'==1);
  [slackrow,~,slackcoef] = find(prob.A(:,slackcol));
  
  change = true;
  while change
    change = false;
    
    % Repair by increasing slack
    pfeas   = prob.A(slackrow,:) *self.xx - prob.b(slackrow)*self.tau;
    stepdir = max(-pfeas ./ slackcoef, 0);
    stepdir(abs(stepdir) <= 1e-5) = 0;
    if any(stepdir)
      change = true;
      self.xx(slackcol) = self.xx(slackcol) + stepdir - sign(stepdir)*1e-6;
    end
    
    % Repair by decreasing slack (don't accept half-way repairs)
    pfeas   = prob.A(slackrow,:) *self.xx - prob.b(slackrow)*self.tau;
    stepdir = min(-pfeas ./ slackcoef, 0);
    stepdir(self.xx(slackcol) + stepdir < 0) = 0;
    stepdir(abs(stepdir) <= 1e-5) = 0;
    if any(stepdir)
      change = true;
      self.xx(slackcol) = self.xx(slackcol) + stepdir - sign(stepdir)*1e-6;
    end
    
    % Repair by unique repair-direction (disregarding slack)
    unirows  = (sum(prob.A(slackrow,:)' ~= 0)'==2);
    unislackcol  = slackcol(unirows);
    unislackrow  = slackrow(unirows);
    unislackcoef = slackcoef(unirows);
    
    notslack = true(size(prob.A,2),1);
    notslack(unislackcol) = false;
    [univarcol,~,univarcoef] = find(prob.A(unislackrow, notslack)');
    
    pfeas   = prob.A(unislackrow,:) *self.xx - prob.b(unislackrow)*self.tau;
    stepdir1 = max(min(-pfeas ./ unislackcoef, 0), -self.xx(unislackcol));
    stepdir1(abs(stepdir1) <= 1e-5) = 0;
    
    pfeas = pfeas + unislackcoef.*stepdir1;
    stepdir2 = -pfeas ./ univarcoef;
    stepdir2(abs(stepdir2) <= 1e-5) = 0;
    
    if any(stepdir1) || any(stepdir1)
      change = true;
      self.xx(unislackcol) = self.xx(unislackcol) + stepdir1 - sign(stepdir1)*1e-6;
      self.xx(univarcol)   = self.xx(univarcol) + stepdir2 - sign(stepdir2)*1e-6;
    end
    
  end

end

