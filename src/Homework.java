import java.util.ArrayList;
import java.util.Random;


public class Homework {
	private static final Random rnd = new Random(222);
	
	// Example
	public static void task0() {
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// construct dummy policy that always executes action 4
		ArrayList<double[][]> policies = new ArrayList<double[][]>();
		double[][] policy = new double[cmdp.getNumStates()][cmdp.getNumActions()];
		for(int s=0; s<cmdp.getNumStates(); s++) {
			policy[s][4] = 1.0;
		}
		policies.add(policy);
		
		// construct dummy solution object with dummy expectations 0.0
		Solution solution = new Solution(policies, 0.0, 0.0, new double[]{0.0}, new double[]{0.0});
		
		// use the simulator to execute on run
		Simulator sim = new Simulator(rnd);
		sim.printActions();
		sim.simulate(cmdps, solution, 1);
	}
	
	// Solve unconstrained problem for 1 agent with value iteration
	public static void task1() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Solve the problem without constraints
		PlanningAlgorithm alg = new PlanningAlgorithm();		
		Solution solution = alg.solveVI(cmdps, 0.95);
		System.out.println("Expected reward: "+solution.getExpectedReward());
		System.out.println("Expected cost: "+solution.getExpectedCost());
		
		// Simulate solution
		System.out.println();
		Simulator sim = new Simulator(rnd);
		sim.simulate(cmdps, solution, 1000);
		
		// Print policy of agent 0
		int agentID = 0;
		double[][] policy = solution.getPolicy(agentID);
		System.out.println();
		for(int s=0; s<cmdps[agentID].getNumStates(); s++) {
			System.out.print("State "+s+": ");
			for(int a=0; a<cmdps[agentID].getNumActions(); a++) {
				System.out.print(policy[s][a]+" ");
			}
			System.out.println();
		}
		
	}
	
	// Solve unconstrained problem for 1 agent with cost
	public static void task2() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Assign cost
		cmdp.assignCost(0, 0, 0); // TODO add costs to the state action pairs
		
		// Solve the problem without constraints
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Solution solution = alg.solveUnconstrained(cmdps);
		System.out.println("Expected reward: "+solution.getExpectedReward());
		System.out.println("Expected cost: "+solution.getExpectedCost());
		
		// Simulate solution
		System.out.println();
		Simulator sim = new Simulator(rnd);
		sim.simulate(cmdps, solution, 1000);
		
		// Print policy of agent 0
		int agentID = 0;
		double[][] policy = solution.getPolicy(agentID);
		System.out.println();
		for(int s=0; s<cmdps[agentID].getNumStates(); s++) {
			System.out.print("State "+s+": ");
			for(int a=0; a<cmdps[agentID].getNumActions(); a++) {
				System.out.print(policy[s][a]+" ");
			}
			System.out.println();
		}
		
	}
	
	// Solve constrained problem for 1 agent
	public static void task3() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Assign cost
		cmdp.assignCost(0, 0, 0); // TODO add costs to the state action pairs
		
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Solution solution = alg.solve(cmdps, 20.0);
		double expectedReward = solution.getExpectedReward();
		System.out.println("Expected reward budget 20: "+expectedReward);
		
		// TODO print expected reward as function of cost limit L
	}
	
	// Solve constrained problem for 2 agents with trivial budget split
	public static void task4() {		
		// Get CMDP models
		CMDP cmdpChild = UserGenerator.getCMDPChild();
		CMDP cmdpAdult = UserGenerator.getCMDPAdult();
		
		// Assign cost to child
		cmdpChild.assignCost(0, 0, 0); // TODO add costs to the state action pairs
		
		// Assign cost to adult
		cmdpAdult.assignCost(0, 0, 0); // TODO add costs to the state action pairs
		
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Simulator sim = new Simulator(rnd);
		
		// Solve both problems separately without constraints and print expectations
		System.out.println("=========== UNCONSTRAINED ===========");
		for(int i=0; i<2; i++) {
			CMDP cmdp = (i==0) ? cmdpChild : cmdpAdult;
			Solution sol = alg.solveUnconstrained(new CMDP[]{cmdp});
			double expectedReward0 = sol.getExpectedReward();
			double expectedCost0 = sol.getExpectedCost();
			System.out.println("Expected reward agent "+i+": "+expectedReward0);
			System.out.println("Expected cost agent "+i+": "+expectedCost0);
		}
		
		// trivial budget split: invest 10 in each agent
		System.out.println();
		System.out.println("=========== SEPARATE PLANNING ===========");
		
		double expectedReward = 0.0;
		double expectedCost = 0.0;
		for(int i=0; i<2; i++) {
			CMDP cmdp = (i==0) ? cmdpChild : cmdpAdult;
			Solution sol = alg.solve(new CMDP[]{cmdp}, 99999.0); // TODO replace the number with the correct limit
			double expectedReward0 = sol.getExpectedReward();
			double expectedCost0 = sol.getExpectedCost();
			System.out.println("Expected reward agent "+i+": "+expectedReward0);
			System.out.println("Expected cost agent "+i+": "+expectedCost0);
			expectedReward += expectedReward0;
			expectedCost += expectedCost0;
		}
		System.out.println("Expected reward: "+expectedReward);
		System.out.println("Expected cost: "+expectedCost);
		
		// multi-agent problem: invest 20 in total
		Solution combinedSolution = alg.solve(new CMDP[]{cmdpChild, cmdpAdult}, 99999.0); // TODO replace the number with the correct limit
		System.out.println();
		System.out.println("=========== MULTI-AGENT PLANNING ===========");
		System.out.println("Expected reward: "+combinedSolution.getExpectedReward());
		System.out.println("Expected reward agent 0: "+combinedSolution.getExpectedReward(0));
		System.out.println("Expected reward agent 1: "+combinedSolution.getExpectedReward(1));
		System.out.println("Expected cost total: "+combinedSolution.getExpectedCost());
		System.out.println("Expected cost agent 0: "+combinedSolution.getExpectedCost(0));
		System.out.println("Expected cost agent 1: "+combinedSolution.getExpectedCost(1));
		
		// simulate
		sim.simulate(new CMDP[]{cmdpChild, cmdpAdult}, combinedSolution, 10000);
	}
	
	public static void main(String[] args) {
		task1();
	}
}
