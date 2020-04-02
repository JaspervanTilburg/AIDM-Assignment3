import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.*;


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
		for (int s = 0; s < cmdp.getNumStates(); s++) {
			for (int a = 0; a < cmdp.getNumActions(); a++) {
				cmdp.assignCost(s, a, 2*a);
			}
		}

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
		for (int s = 0; s < cmdp.getNumStates(); s++) {
			for (int a = 0; a < cmdp.getNumActions(); a++) {
				cmdp.assignCost(s, a, 2*a);
			}
		}

		PlanningAlgorithm alg = new PlanningAlgorithm();
		for (int budget = 0; budget <= 50; budget += 5) {
            Solution solution = alg.solve(cmdps, budget);
            double expectedReward = solution.getExpectedReward();
            System.out.println("Expected reward budget " + budget + ": " + expectedReward);
        }
	}
	
	// Solve constrained problem for 2 agents with trivial budget split
	public static void task4() {		
		// Get CMDP models
		CMDP cmdpChild = UserGenerator.getCMDPChild();
		CMDP cmdpAdult = UserGenerator.getCMDPAdult();
		
		// Assign cost to child
		for (int s = 0; s < cmdpChild.getNumStates(); s++) {
			for (int a = 0; a < cmdpChild.getNumActions(); a++) {
				cmdpChild.assignCost(s, a, 2*a);
			}
		}
		
		// Assign cost to adult
		for (int s = 0; s < cmdpAdult.getNumStates(); s++) {
			for (int a = 0; a < cmdpAdult.getNumActions(); a++) {
				cmdpAdult.assignCost(s, a, 2*a);
			}
		}
		
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Simulator sim = new Simulator(rnd);
		double budget = 20;
		
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
			Solution sol = alg.solve(new CMDP[]{cmdp}, budget/2);
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
		Solution combinedSolution = alg.solve(new CMDP[]{cmdpChild, cmdpAdult}, budget);
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

	public static void task5() {
		// Only children
		System.out.println("Experiment with varying number of children");
		for (int i = 2; i <= 50; i++) {
			runExperiment(i, 0, i*10);
		}

		// Only adults
		System.out.println("Experiment with varying number of adults");
		for (int i = 2; i <= 50; i++) {
			runExperiment(0, i, i*10);
		}

		// Children and adults equally divided.
		System.out.println("Experiment with equally divided number of children and adults");
		for (int i = 1; i <= 25; i++) {
			runExperiment(i, i, i*10);
		}

		// All compositions of 20 agents
		System.out.println("Experiment with varying composition of children and adults");
		for (int i = 0; i <= 50; i++) {
			runExperiment(i, 50 - i, 200);
		}
	}

	private static void runExperiment(int numChildAgents, int numAdultAgents, double L) {
		CMDP[] agents = new CMDP[numChildAgents + numAdultAgents];

		// Generate child users and assign costs
		for (int i = 0; i < numChildAgents; i++) {
			CMDP cmdpChild = UserGenerator.getCMDPChild();
			agents[i] = cmdpChild;
			for (int s = 0; s < cmdpChild.getNumStates(); s++) {
				for (int a = 0; a < cmdpChild.getNumActions(); a++) {
					cmdpChild.assignCost(s, a, 2*a);
				}
			}
		}

		// Generate adult users and assign costs
		for (int i = 0; i < numAdultAgents; i++) {
			CMDP cmdpAdult = UserGenerator.getCMDPAdult();
			agents[numChildAgents + i] = cmdpAdult;
			for (int s = 0; s < cmdpAdult.getNumStates(); s++) {
				for (int a = 0; a < cmdpAdult.getNumActions(); a++) {
					cmdpAdult.assignCost(s, a, 2*a);
				}
			}
		}

		// Setup the solver as a timed callable task
		PlanningAlgorithm alg = new PlanningAlgorithm();
		ExecutorService executor = Executors.newCachedThreadPool();
		Callable<Object> task = () -> alg.solve(agents, L);
		Future<Object> future = executor.submit(task);

		// Try to run the solver task with a timeout of 2 minutes
		try {
			long start = System.currentTimeMillis();
			Solution sol = (Solution) future.get(2, TimeUnit.MINUTES);
			long end = System.currentTimeMillis();
			System.out.println("Expected reward: " + sol.getExpectedReward());
			System.out.println("Expected cost: " + sol.getExpectedCost());
			System.out.println("Runtime in ms: " + (end - start));
		} catch (TimeoutException ex) {
			System.out.println("Timeout");
		} catch (InterruptedException e) {
			System.out.println("InterruptedException");
		} catch (ExecutionException e) {
			System.out.println("ExecutionException");
		}
	}
	
	public static void main(String[] args) {
//		task0();
//		task1();
//		task2();
//		task3();
//		task4();
		task5();
	}
}
