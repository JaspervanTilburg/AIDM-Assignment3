import java.util.ArrayList;


public class Solution {
	private ArrayList<double[][]> policies;
	private double expectedReward;
	private double expectedCost;
	private double[] expectedRewardAgent;
	private double[] expectedCostAgent;
	
	public Solution(ArrayList<double[][]> policies, double expectedReward, double expectedCost, double[] expectedRewardAgent, double[] expectedCostAgent) {
		this.policies = policies;
		this.expectedReward = expectedReward;
		this.expectedCost = expectedCost;
		this.expectedRewardAgent = expectedRewardAgent;
		this.expectedCostAgent = expectedCostAgent;
	}
	
	public double[][] getPolicy(int i) {
		return policies.get(i);
	}
	
	public double getExpectedReward() {
		return expectedReward;
	}
	
	public double getExpectedCost() {
		return expectedCost;
	}
	
	public double getExpectedReward(int i) {
		return expectedRewardAgent[i];
	}
	
	public double getExpectedCost(int i) {
		return expectedCostAgent[i];
	}
	
}
