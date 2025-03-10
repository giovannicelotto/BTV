# %%
import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_geometric.datasets import KarateClub
from torch_geometric.data import DataLoader
# %%
# Load the Karate Club dataset
dataset = KarateClub()
data = dataset[0]  # The dataset has only one graph
data.test_mask = ~data.train_mask  # Inverse of train_mask
# %%
# Define a simple Graph Neural Network
class GNN(torch.nn.Module):
    def __init__(self):
        super(GNN, self).__init__()
        self.conv1 = GCNConv(dataset.num_features, 16)  # First graph convolution
        self.conv2 = GCNConv(16, dataset.num_classes)   # Second graph convolution
    
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1)
# %%
# Initialize model, optimizer, and loss function
model = GNN()
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
loss_fn = torch.nn.NLLLoss()

def train():
    model.train()
    optimizer.zero_grad()
    out = model(data)
    loss = loss_fn(out[data.train_mask], data.y[data.train_mask])
    loss.backward()
    optimizer.step()
    return loss.item()

# Training loop
for epoch in range(200):
    loss = train()
    if epoch % 20 == 0:
        print(f"Epoch {epoch}, Loss: {loss:.4f}")
# %%
# Evaluate the model
model.eval()
pred = model(data).argmax(dim=1)
correct = (pred[data.test_mask] == data.y[data.test_mask]).sum()
accuracy = int(correct) / int(data.test_mask.sum())
print(f"Test Accuracy: {accuracy:.4f}")

# %%
import torch.nn as nn
import torch.nn.functional as F
class DNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(DNN, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)  # First hidden layer
        self.fc2 = nn.Linear(hidden_dim, output_dim)  # Output layer

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return F.log_softmax(x, dim=1)
input_dim = dataset.num_features  # Feature dimension
hidden_dim = 16  # Hidden layer size
output_dim = dataset.num_classes  # Number of classes

# Initialize model, optimizer, and loss function
model = DNN(input_dim, hidden_dim, output_dim)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
loss_fn = nn.NLLLoss()

# Define a test mask (inverse of train_mask)
data.test_mask = ~data.train_mask  # Nodes not used for training

# Training function
def train():
    model.train()
    optimizer.zero_grad()
    out = model(data.x)  # Use only node features, no graph structure
    loss = loss_fn(out[data.train_mask], data.y[data.train_mask])
    loss.backward()
    optimizer.step()
    return loss.item()

# Training loop
for epoch in range(300):
    loss = train()
    if epoch % 20 == 0:
        print(f"Epoch {epoch}, Loss: {loss:.4f}")

# Evaluation
model.eval()
pred = model(data.x).argmax(dim=1)  # Predict class labels
correct = (pred[data.test_mask] == data.y[data.test_mask]).sum()
accuracy = int(correct) / int(data.test_mask.sum())
print(f"DNN Test Accuracy: {accuracy:.4f}")
# %%
