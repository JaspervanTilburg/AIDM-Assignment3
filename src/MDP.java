
public class MDP {
	private int nStates;
	private int nActions;
	private int initialState;
	
	private double discountFactor;
	
	private double[][] rewardFunction;
	private double minReward = Double.POSITIVE_INFINITY;
	private double maxReward = Double.NEGATIVE_INFINITY;
	
	private double[][][] transitionFunction;
	
	
	public MDP(int nStates, int nActions, int initialState, double discountFactor) {
		this.nStates = nStates;
		this.nActions = nActions;
		this.initialState = initialState;
		this.discountFactor = discountFactor;
	}
	
	/**
	 * Get number of states
	 * @return number of states
	 */
	public int getNumStates() {
		return nStates;
	}
	
	/**
	 * Get number of actions
	 * @return number of actions
	 */
	public int getNumActions() {
		return nActions;
	}
	
	/**
	 * Get discount factor
	 * @return discount factor
	 */
	public double getDiscountFactor() {
		return discountFactor;
	}
	
	/**
	 * Get initial state
	 * @return initial state
	 */
	public int getInitialState() {
		return initialState;
	}
	
	/**
	 * Set reward function
	 * @param rewardFunction reward function
	 */
	public void setRewardFunction(double[][] rewardFunction) {
		this.rewardFunction = rewardFunction;
		
		minReward = Double.POSITIVE_INFINITY;
		maxReward = Double.NEGATIVE_INFINITY;
		
		for(int s=0; s<nStates; s++) {
			for(int a=0; a<nActions; a++) {
				minReward = Math.min(minReward, rewardFunction[s][a]);
				maxReward = Math.max(maxReward, rewardFunction[s][a]);
			}
		}
	}
	
	/**
	 * Get reward function
	 * @return reward function
	 */
	public double[][] getRewardFunction() {
		return rewardFunction;
	}
	
	/**
	 * Get reward R(s,a)
	 * @param s state s
	 * @param a action a
	 * @return reward R(s,a)
	 */
	public double getReward(int s, int a) {
		assert s<nStates && a<nActions;
		return rewardFunction[s][a];
	}
	
	/**
	 * Get minimum instantaneous reward
	 * @return min reward
	 */
	public double getMinReward() {
		return minReward;
	}
	
	/**
	 * Get maximum instantaneous reward
	 * @return max reward
	 */
	public double getMaxReward() {
		return maxReward;
	}
	
	/**
	 * Set transition function
	 * @param transitionFunction transition function
	 */
	public void setTransitionFunction(double[][][] transitionFunction) {
		this.transitionFunction = transitionFunction;
	}
	
	/**
	 * Get transition probability
	 * @param s state s
	 * @param a action a
	 * @param sNext successor state sNext
	 * @return
	 */
	public double getTransitionProbability(int s, int a, int sNext) {
		return transitionFunction[s][a][sNext];
	}
}
