
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class ProbabilitySample {
	private List<Item> items;
	private Random rnd;
	private double probabilitySum = 0.0;
	
	public ProbabilitySample(Random rnd) {
		items = new ArrayList<Item>();
		this.rnd = rnd;
	}
	
	public void addItem(int item, double probability) {
		assert probability >= 0.0 && probability <= 1.0 : "PROB: "+probability;		
		if(probability > 0.0) {
			items.add(new Item(item,probability));
			probabilitySum += probability;
		}
	}
	
	public int sampleItem() {
		assert Math.abs(probabilitySum-1.0) < 0.001 : "No valid probability distribution: "+probabilitySum;
		assert items.size() > 0 : "No items added";
		
		double cumulative = 0.0;
		double randomNumber = rnd.nextDouble();
		int retItem = items.get(items.size()-1).item;
		
		for(Item item : items) {
			cumulative += item.probability;
			
			if(randomNumber <= cumulative) {
				retItem = item.item;
				break;
			}
		}
		
		return retItem;
	}
	
	private class Item {
		public int item;
		public double probability;
		
		public Item(int item, double probability) {
			this.item = item;
			this.probability = probability;
		}
	}
}
