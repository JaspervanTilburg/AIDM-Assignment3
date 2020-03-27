

public class CMDP extends MDP {	
	private double[][] costFunction;
	
	public CMDP(int nStates, int nActions, int initialState, double discountFactor, double[][] costFunction) {
		super(nStates, nActions, initialState, discountFactor);
		this.costFunction = costFunction;
	}
	
	/**
	 * Get cost / resource consumption for resource k when executing action a in state s
	 * @param k resource id
	 * @param s state
	 * @param a action
	 * @return cost / resource consumption
	 */
	public double getCost(int s, int a) {
		assert s<super.getNumStates() && a<super.getNumActions();
		return costFunction[s][a];
	}
	
	/**
	 * Get cost function
	 * @return cost function
	 */
	public double[][] getCostFunction() {
		return costFunction;
	}
	
	/**
	 * Assign cost to a state-action pair
	 * @param s state
	 * @param a action
	 * @param cost cost corresponding to s,a
	 */
	public void assignCost(int s, int a, double cost) {
		costFunction[s][a] = cost;
	}
}
