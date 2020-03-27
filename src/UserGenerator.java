import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class UserGenerator {
	
	public static CMDP getCMDPChild() {
		CMDP cmdp = null;
		int[][][] Tstate;
		double[][][] Tprob;
		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader("synthetic_ad.txt"));

			int initialState = 0;
			int numStates = Integer.valueOf(reader.readLine());
			int numActions = Integer.valueOf(reader.readLine());

			double[][] costFunction = new double[numStates][numActions];
			Tstate = new int[numStates][numActions][];
			Tprob = new double[numStates][numActions][];
			
			double[][] rewardFunction = new double[numStates][numActions];
			double[][][] transitionFunction = new double[numStates][numActions][numStates];

			// discount factor, unused
			reader.readLine();

			// repeated action description.
			for (int a = 0; a < numActions; a++) {
				int actionID = Integer.valueOf(reader.readLine());
				
				if (a != actionID) {
					throw new RuntimeException("Unexpected action ID (" + a + " != " + actionID + ").");
				}
				
				// define reward
				rewardFunction[9][a] = 200.0;

				String[] chunks = null;
				for (int s = 0; s < numStates; s++) {
					chunks = reader.readLine().split(" ");

					int stateID = Integer.valueOf(chunks[0]);
					if (s != stateID) {
						throw new RuntimeException("Unexpected state ID (" + s + " != " + stateID + ").");
					}

					int numDests = (chunks.length-1)/2;
					Tstate[s][a] = new int[numDests];
					Tprob[s][a] = new double[numDests];

					for (int i = 0; i < numDests; i++) {
						int sNext = Integer.valueOf(chunks[1+(i*2)].substring(1));
						double prob = Double.valueOf(chunks[2+(i*2)].substring(0, chunks[2+(i*2)].length()-1));

						Tstate[s][a][i] = sNext;
						Tprob[s][a][i] = prob;
						
						transitionFunction[s][a][sNext] = prob;
					}
				}

				// Rewards, unused.
				reader.readLine();

				// Resource costs, unused.
				reader.readLine();
				
			}
			
			cmdp = new CMDP(numStates, numActions, initialState, 0.95, costFunction);
			cmdp.setRewardFunction(rewardFunction);
			cmdp.setTransitionFunction(transitionFunction);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) { }
			}
		}
		
		return cmdp;
	}
	
	public static CMDP getCMDPAdult() {
		CMDP cmdp = null;
		int[][][] Tstate;
		double[][][] Tprob;
		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader("synthetic_ad.txt"));

			int initialState = 0;
			int numStates = Integer.valueOf(reader.readLine());
			int numActions = Integer.valueOf(reader.readLine());

			double[][] costFunction = new double[numStates][numActions];
			Tstate = new int[numStates][numActions][];
			Tprob = new double[numStates][numActions][];
			
			double[][] rewardFunction = new double[numStates][numActions];
			double[][][] transitionFunction = new double[numStates][numActions][numStates];

			// discount factor, unused
			reader.readLine();

			// repeated action description.
			for (int a = 0; a < numActions; a++) {
				int actionID = Integer.valueOf(reader.readLine());
				
				if (a != actionID) {
					throw new RuntimeException("Unexpected action ID (" + a + " != " + actionID + ").");
				}
				
				// define reward
				rewardFunction[9][a] = 1000.0;

				String[] chunks = null;
				for (int s = 0; s < numStates; s++) {
					chunks = reader.readLine().split(" ");

					int stateID = Integer.valueOf(chunks[0]);
					if (s != stateID) {
						throw new RuntimeException("Unexpected state ID (" + s + " != " + stateID + ").");
					}

					int numDests = (chunks.length-1)/2;
					Tstate[s][a] = new int[numDests];
					Tprob[s][a] = new double[numDests];

					for (int i = 0; i < numDests; i++) {
						int sNext = Integer.valueOf(chunks[1+(i*2)].substring(1));
						double prob = Double.valueOf(chunks[2+(i*2)].substring(0, chunks[2+(i*2)].length()-1));

						Tstate[s][a][i] = sNext;
						Tprob[s][a][i] = prob;
						
						transitionFunction[s][a][sNext] = prob;
					}
				}

				// Rewards, unused.
				reader.readLine();

				// Resource costs, unused.
				reader.readLine();
				
			}
			
			// change transitions for state 7
			for(int a=0; a<numActions; a++) {
				double increment = 0.7;
				transitionFunction[7][a][7] += increment;
				
				for(int sNext=0; sNext<numStates; sNext++) {
					transitionFunction[7][a][sNext] = transitionFunction[7][a][sNext] / (1.0+increment);
				}
			}
			
			cmdp = new CMDP(numStates, numActions, initialState, 0.95, costFunction);
			cmdp.setRewardFunction(rewardFunction);
			cmdp.setTransitionFunction(transitionFunction);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) { }
			}
		}
		
		return cmdp;
	}
}
